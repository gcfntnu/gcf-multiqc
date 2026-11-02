"""Installed reporting commands used by the suite; no BFQ or workflow execution.

The invocation and sequencing config shapes follow gcf-bfq 2b1fbc5 and
gcf-workflows 26fafe4, with BFQ issue #173's deduplicated multi-pass selection.
All data below are synthetic.
"""

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml


FIXTURES = Path(__file__).parent / "fixtures"


def run_report(command, cwd):
    result = subprocess.run(
        [sys.executable, "-m", "multiqc", "--strict", "--no-version-check", "--no-ai", *command],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def test_workflow_relative_filename_keeps_report_and_exports_under_interim_directory(tmp_path, general_stats_cells):
    interim = Path("data/tmp/singlecell/bfq")
    logs = tmp_path / interim / "logs"
    figures = tmp_path / interim / "figs"
    shutil.copytree(FIXTURES / "compatibility", logs)
    shutil.copytree(FIXTURES / "singlecell/parse", logs / "parse")
    figures.mkdir()
    (figures / "synthetic_mqc.yaml").write_text(
        yaml.safe_dump(
            {
                "id": "synthetic_summary",
                "plot_type": "html",
                "section_name": "Synthetic summary",
                "data": "<p>Summary available.</p>",
            }
        )
    )
    config = yaml.safe_load((FIXTURES / "profiles/multiqc_config-singlecell.yaml").read_text())
    config["custom_logo"] = str(FIXTURES / "profiles/logo_ntnu.png")
    config["custom_data"] = {
        "general_statistics": {
            "plot_type": "generalstats",
            "pconfig": [{"External_ID": {"namespace": "gcf", "format": "{}"}}],
            "data": {"001": {"External_ID": "0007"}},
        }
    }
    config_path = interim / ".multiqc_config.yaml"
    (tmp_path / config_path).write_text(yaml.safe_dump(config))
    report_path = interim / "multiqc_SYNTHETIC.html"
    modules = [
        "cellranger",
        "starsolo",
        "fastp",
        "fastqc_rnaseq",
        "custom_content",
        "cellranger_count",
        "kraken",
        "parse",
    ]
    command = ["--config", str(config_path)]
    for module in modules:
        command += ["-m", module]
    command += [
        "--filename",
        str(report_path),
        "-f",
        "--interactive",
        "--no-megaqc-upload",
        str(interim / "logs"),
        str(interim / "figs"),
    ]
    run_report(command, tmp_path)
    assert (tmp_path / report_path).stat().st_size > 1000
    assert not (tmp_path / "multiqc_report.html").exists()
    data = json.loads((tmp_path / interim / "multiqc_SYNTHETIC_data/multiqc_data.json").read_text())
    assert data["report_saved_raw_data"]["multiqc_parse"]["001"]["number_of_reads"] == 1000
    assert data["report_general_stats_data"]["custom_content"]["001"]["External_ID"] == "0007"
    assert data["report_saved_raw_data"]["multiqc_general_stats"]["001"]["custom_content_gcf-External_ID"] == "0007"
    assert general_stats_cells((tmp_path / report_path).read_text(), "custom_content_gcf-External_ID")["001"] == "0007"
    assert "Synthetic summary" in (tmp_path / report_path).read_text()


@pytest.mark.parametrize("demultiplexer", ["bclconvert", "bcl2fastq", "custom_content"])
def test_bfq_sequencing_invocation_keeps_custom_missing_and_zero_assignments(tmp_path, demultiplexer):
    # Planned separate conversion passes select custom_content for the combined report.
    report_dir = tmp_path / "sequencing report"
    inputs = report_dir / "inputs"
    inputs.mkdir(parents=True)
    for module in ["interop", *([] if demultiplexer == "custom_content" else [demultiplexer])]:
        for source in (FIXTURES / "sequencing" / module).iterdir():
            shutil.copy(source, inputs / source.name)
    (inputs / "bfq_sequencing_summary_mqc.yaml").write_text(
        yaml.safe_dump(
            {
                "id": "bfq_sequencing_summary",
                "section_name": "BFQ sequencing summary",
                "plot_type": "html",
                "data": "<p>Synthetic run summary.</p>",
            }
        )
    )
    assignments = {
        "SYNTHETIC / 001": {"project": "SYNTHETIC", "sample_id": "001", "expected": "Yes", "assigned_reads": 0},
        "SYNTHETIC / NA": {
            "project": "SYNTHETIC",
            "sample_id": "NA",
            "expected": "Yes",
            "assigned_reads": "Unavailable",
        },
    }
    config = {
        "title": "Sequencing QC — SYNTHETIC",
        "megaqc_upload": False,
        "megaqc_url": None,
        "megaqc_access_token": None,
        "custom_data": {
            "bfq_sample_assignment": {
                "section_name": "Expected samples and demultiplexing assignment",
                "plot_type": "table",
                "description": "Missing statistics are unavailable, not zero. Explicit zero counts are retained.",
                "pconfig": {
                    "id": "bfq_sample_assignment",
                    "title": "Sample assignment",
                    "col1_header": "Project / sample",
                },
                "headers": {
                    "project": {"title": "Project"},
                    "sample_id": {"title": "Sample ID"},
                    "expected": {"title": "In SampleSheet"},
                    "assigned_reads": {"title": "Assigned clusters", "format": "{:,.0f}"},
                },
                "data": assignments,
            }
        },
        "module_order": list(dict.fromkeys(["custom_content", "interop", demultiplexer])),
    }
    config_path = report_dir / "multiqc_config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    filename = "sequencer_stats_SYNTHETIC.html"
    command = [
        "-q",
        "--interactive",
        "--force",
        "--no-megaqc-upload",
        "--config",
        str(config_path),
        "--outdir",
        str(report_dir),
        "--filename",
        filename,
    ]
    for module in config["module_order"]:
        command += ["-m", module]
    run_report([*command, str(inputs)], report_dir)
    html = (report_dir / filename).read_text()
    assert "Synthetic run summary." in html
    report = json.loads((report_dir / "sequencer_stats_SYNTHETIC_data/multiqc_data.json").read_text())
    raw = report["report_saved_raw_data"]
    assert raw["multiqc_bfq_sample_assignment"] == assignments
    assert set(key for key in raw if key.startswith("multiqc_bfq_sample_assignment")) == {
        "multiqc_bfq_sample_assignment"
    }
    assert "interop_runsummary" in raw
    if demultiplexer != "custom_content":
        assert f"multiqc_{demultiplexer}_bysample" in raw
