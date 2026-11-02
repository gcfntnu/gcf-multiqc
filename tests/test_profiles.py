"""Exercise the existing seven report selections without running any workflow."""

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

FIXTURES = Path(__file__).parent / "fixtures"
PROFILES = {
    "rnaseq": "fastq_screen,star,picard,fastp,fastqc_rnaseq,salmon,rseqc,custom_content,kraken",
    "microbiome": "fastq_screen,fastp,fastqc_rnaseq,custom_content,qiime2",
    "singlecell": "cellranger,starsolo,fastp,fastqc_rnaseq,custom_content,cellranger_count,kraken,parse",
    "smallrna": "star,picard,fastp,fastqc_rnaseq,unitas,custom_content",
    "default": "fastq_screen,fastp,fastqc_rnaseq,custom_content,kraken",
    "wgs": "fastq_screen,fastqc_rnaseq,qualimap,fastp,picard,kraken",
    "metagenome": "fastqc_rnaseq,fastp,kraken,bracken,custom_content",
}


@pytest.mark.parametrize("profile", PROFILES)
def test_existing_workflow_profile_accepts_legacy_selection(tmp_path, profile, general_stats_cells):
    inputs = tmp_path / "inputs"
    shutil.copytree(FIXTURES / "compatibility", inputs)
    custom_fixtures = {
        "microbiome": "microbiome_smallrna/qiime2",
        "smallrna": "microbiome_smallrna/unitas",
        "singlecell": "singlecell",
        "metagenome": "bracken",
    }
    if profile in custom_fixtures:
        shutil.copytree(FIXTURES / custom_fixtures[profile], inputs / "specific")
    config = yaml.safe_load((FIXTURES / "profiles" / f"multiqc_config-{profile}.yaml").read_text())
    config.update(
        custom_logo=str(FIXTURES / "profiles/logo_ntnu.png"),
        no_version_check=True,
        megaqc_upload=False,
        ai_summary=False,
        ai_summary_full=False,
    )
    config["custom_data"] = {
        "general_statistics": {
            "plot_type": "generalstats",
            "pconfig": [{"External_ID": {"namespace": "gcf", "format": "{}"}}],
            "data": {"Synthetic": {"External_ID": "0007"}},
        }
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    output = tmp_path / "output"
    cmd = [
        sys.executable,
        "-m",
        "multiqc",
        "--strict",
        "--interactive",
        "--no-version-check",
        "--no-ai",
        "--no-megaqc-upload",
        "-c",
        str(config_path),
        "-o",
        str(output),
    ]
    for module in PROFILES[profile].split(","):
        cmd += ["-m", module]
    cmd.append(str(inputs))
    result = subprocess.run(
        cmd, cwd=tmp_path, capture_output=True, text=True, timeout=60, env={**os.environ, "TMPDIR": str(tmp_path)}
    )
    assert result.returncode == 0, result.stdout + result.stderr
    data = json.loads((output / "multiqc_data/multiqc_data.json").read_text())
    assert data["report_general_stats_data"]["fastp"]["Synthetic"]["before_filtering_total_reads"] == 1000
    assert "fastqc" not in data["report_general_stats_data"]
    assert "fastqc_rnaseq" not in data["report_general_stats_data"]
    if "custom_content" in PROFILES[profile].split(","):
        assert data["report_general_stats_data"]["custom_content"]["Synthetic"]["External_ID"] == "0007"
        assert (
            data["report_saved_raw_data"]["multiqc_general_stats"]["Synthetic"]["custom_content_gcf-External_ID"]
            == "0007"
        )
        assert (
            general_stats_cells((output / "multiqc_report.html").read_text(), "custom_content_gcf-External_ID")[
                "Synthetic"
            ]
            == "0007"
        )
        with (output / "multiqc_data/multiqc_general_stats.txt").open() as handle:
            rows = {row["Sample"]: row for row in csv.DictReader(handle, delimiter="\t")}
        assert rows["Synthetic"]["custom_content_gcf-External_ID"] == "0007"
    raw = data["report_saved_raw_data"]
    if profile == "microbiome":
        assert raw["multiqc_qiime2_dada2"]["Total"]["001"]["passed"] == 105
        assert raw["multiqc_qiime2_rpca"]["001"]["name"] == "001"
    elif profile == "smallrna":
        assert raw["multiqc_unitas"]["001"]["miRNA Percentage"] == 0.1
        assert raw["multiqc_unitas"]["1"]["miRNA Percentage"] == 0.1
    elif profile == "singlecell":
        assert raw["multiqc_parse"]["all_samples"]["number_of_reads"] == 2000
        assert raw["multiqc_starsolo"]["001"]["Mean Reads per Cell"] == 85.5
        assert raw["multiqc_cellranger_count"]["001"]["Number of Reads"] == 1000
        assert raw["multiqc_cellranger"]["001"]["number_reads"] == 1000
    elif profile == "metagenome":
        assert raw["multiqc_bracken"] and raw["multiqc_kraken"]
        assert raw["multiqc_bracken"] != raw["multiqc_kraken"]
    assert (output / "multiqc_report.html").stat().st_size > 1000
