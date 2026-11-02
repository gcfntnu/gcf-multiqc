"""Render public synthetic inputs through installed MultiQC and GCF entry points."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import pytest
import yaml


FIXTURES = Path(__file__).parent / "fixtures" / "compatibility"


def render(tmp_path, inputs, modules, extra_config=None):
    config = {
        "no_version_check": True,
        "megaqc_upload": False,
        "ai_summary": False,
        "ai_summary_full": False,
        "data_dump_file": True,
        "preserve_module_raw_data": True,
        "read_count_multiplier": 1,
        "read_count_prefix": "",
        "read_count_desc": "",
        **(extra_config or {}),
    }
    config_file = tmp_path / "multiqc_test.yaml"
    config_file.write_text(yaml.safe_dump(config))
    output = tmp_path / "output"
    command = [
        sys.executable,
        "-m",
        "multiqc",
        str(inputs),
        "--strict",
        "--require-logs",
        "--no-megaqc-upload",
        "-c",
        str(config_file),
        "-o",
        str(output),
    ]
    for module in modules:
        command.extend(["-m", module])
    completed = subprocess.run(
        command,
        cwd=tmp_path,
        env={**os.environ, "TMPDIR": str(tmp_path)},
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    return (
        json.loads((output / "multiqc_data" / "multiqc_data.json").read_text()),
        (output / "multiqc_report.html").read_text(),
        output,
    )


def rendered_columns(data, section):
    """The rendered table applies placement/visibility after raw header export."""
    return data["report_plot_data"]["general_stats_table"]["datasets"][0]["dt"]["section_by_id"][section][
        "column_by_key"
    ]


@pytest.mark.parametrize("zipped", [False, True])
def test_fastqc_alias_preserves_plots_and_removal_selectors(tmp_path, zipped):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    report_file = FIXTURES / "fastqc" / "fastqc_data.txt"
    if zipped:
        with zipfile.ZipFile(inputs / "Synthetic_fastqc.zip", "w") as archive:
            archive.writestr("Synthetic_fastqc/", "")
            archive.write(report_file, "Synthetic_fastqc/fastqc_data.txt")
    else:
        shutil.copy(report_file, inputs / "fastqc_data.txt")
    shutil.copy(FIXTURES / "Synthetic.fastp.json", inputs)
    data, html, output = render(
        tmp_path,
        inputs,
        ["fastqc_rnaseq", "fastp"],
        {"remove_sections": ["fastqc_rnaseq_sequence_counts", "fastqc_rnaseq_status_checks"]},
    )
    assert 'id="fastqc_rnaseq"' in html
    assert 'id="fastqc_rnaseq_per_base_sequence_content"' in html
    assert 'id="fastqc_rnaseq_sequence_counts"' not in html
    assert 'id="fastqc_rnaseq_status_checks"' not in html
    assert set(data["report_general_stats_data"]) == {"fastp"}
    assert data["report_saved_raw_data"]["multiqc_fastqc"]["Synthetic"]["Total Sequences"] == 1000
    assert (output / "multiqc_data" / "multiqc_fastqc.txt").is_file()
    assert data["report_plot_data"]


@pytest.mark.parametrize("selector", ["fastp", "gcf_fastp"])
def test_fastp_preserves_before_filter_counts_and_configured_columns(tmp_path, selector):
    data, html, _ = render(
        tmp_path,
        FIXTURES / "Synthetic.fastp.json",
        [selector],
        {"table_columns_placement": {"fastp": {"before_filtering_total_reads": 110, "pct_surviving": 112}}},
    )
    values = data["report_general_stats_data"]["fastp"]["Synthetic"]
    headers = data["report_general_stats_headers"]["fastp"]
    assert values["before_filtering_total_reads"] == 1000
    assert values["pct_surviving"] == 80
    assert "filtering_result_passed_filter_reads" not in headers
    assert "pct_surviving_reads" not in headers
    columns = rendered_columns(data, "fastp")
    assert columns["before_filtering_total_reads"]["placement"] == 110
    assert columns["pct_surviving"]["placement"] == 112
    raw = data["report_saved_raw_data"]["multiqc_fastp"]["Synthetic"]
    assert raw["summary"]["after_filtering"]["total_reads"] == 800
    assert 'id="fastp-filtered-reads-chart"' in html


def test_fastp_zero_input_has_no_invented_passing_percentage(tmp_path):
    payload = json.loads((FIXTURES / "Synthetic.fastp.json").read_text())
    for values in payload["summary"].values():
        if isinstance(values, dict):
            for key in values:
                values[key] = 0
    for section in ["filtering_result", "adapter_cutting"]:
        payload[section] = {key: 0 for key in payload[section]}
    fixture = tmp_path / "Zero.fastp.json"
    fixture.write_text(json.dumps(payload))
    data, _, _ = render(tmp_path, fixture, ["fastp"])
    values = data["report_general_stats_data"]["fastp"]["Synthetic"]
    assert values["before_filtering_total_reads"] == 0
    assert "pct_surviving" not in values


def star_inputs(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    text = (FIXTURES / "MultiSplice.Log.final.out").read_text()
    (inputs / "MultiSplice.Log.final.out").write_text(text)
    no_splice = text.replace("Total | 120", "Total | 0").replace("(sjdb) | 90", "(sjdb) | 0")
    for old in ["GT/AG | 100", "GC/AG | 10", "AT/AC | 5", "Non-canonical | 5"]:
        no_splice = no_splice.replace(old, old.rsplit(" ", 1)[0] + " 0")
    (inputs / "NoSplice.Log.final.out").write_text(no_splice)
    no_unique = no_splice.replace("reads number | 50", "reads number | 0")
    no_unique = no_unique.replace("reads % | 50.00%", "reads % | 0.00%")
    no_unique = no_unique.replace("too short | 20.00%", "too short | 70.00%")
    (inputs / "NoUnique.Log.final.out").write_text(no_unique)
    return inputs


@pytest.mark.parametrize("selector", ["star", "gcf_star"])
def test_star_preserves_keys_with_read_and_event_semantics(tmp_path, selector):
    data, html, _ = render(tmp_path, star_inputs(tmp_path), [selector])
    stats = data["report_general_stats_data"]["star"]
    headers = data["report_general_stats_headers"]["star"]
    raw = data["report_saved_raw_data"]["multiqc_star"]
    values = stats["MultiSplice"]
    # Too-many-loci reads are excluded, matching upstream, so this is 70 not 80.
    assert values["total_mapped"] == 70
    assert values["total_mapped_percent"] == 70
    assert values["total_spliced_percent"] == 240
    assert values["novel_spliced_percent"] == 25
    assert stats["NoSplice"]["total_spliced_percent"] == 0
    assert stats["NoSplice"]["novel_spliced_percent"] is None
    assert stats["NoUnique"]["total_spliced_percent"] is None
    assert raw["MultiSplice"]["total_mapped"] == 70
    assert raw["NoUnique"]["total_spliced_percent"] is None
    assert "events" in headers["total_spliced_percent"]["title"]
    assert "max" not in headers["total_spliced_percent"] or headers["total_spliced_percent"]["max"] != 100
    assert "other_uniquely_mapped" not in raw["MultiSplice"]
    assert 'id="star_spliced_reads"' in html
    assert data["report_plot_data"]["star_splice_plot"]["pconfig"]["title"] == "STAR: Splice Events"


def test_star_legacy_smallrna_visibility_and_section_removal(tmp_path):
    data, html, _ = render(
        tmp_path,
        star_inputs(tmp_path),
        ["star"],
        {
            "table_columns_visible": {"STAR": {"total_spliced_percent": False, "novel_spliced_percent": False}},
            "remove_sections": ["star_spliced_reads"],
        },
    )
    headers = rendered_columns(data, "star")
    assert headers["total_spliced_percent"]["hidden"] is True
    assert headers["novel_spliced_percent"]["hidden"] is True
    assert 'id="star_spliced_reads"' not in html
