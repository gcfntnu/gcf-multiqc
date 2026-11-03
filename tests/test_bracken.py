"""GCF Kraken and Bracken coexistence using synthetic public report formats."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml


FIXTURES = Path(__file__).parent / "fixtures" / "bracken"


def render(tmp_path, modules, extra_config=None):
    config_file = tmp_path / "multiqc_test.yaml"
    config_file.write_text(
        yaml.safe_dump(
            {
                "no_version_check": True,
                "megaqc_upload": False,
                "ai_summary": False,
                "ai_summary_full": False,
                "data_dump_file": True,
                "preserve_module_raw_data": True,
                # Exact relevant settings from the GCF metagenome template.
                "extra_fn_clean_exts": [".kraken.kreport", ".bracken"],
                "fn_ignore_files": ["*.bracken_kreport"],
                "kraken": {"top_n": 10},
                "remove_sections": ["kraken-top-n"],
                "table_columns_visible": {"kraken": {"pct_top_one": False, "pct_top_n": False}},
                "table_columns_placement": {
                    "bracken": {"pct_top_one": 200, "pct_top_n": 210},
                    "kraken": {"pct_unclassified": 220},
                },
                **(extra_config or {}),
            }
        )
    )
    output = tmp_path / "output"
    command = [
        sys.executable,
        "-m",
        "multiqc",
        str(tmp_path / "inputs"),
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
    data = json.loads((output / "multiqc_data" / "multiqc_data.json").read_text())
    return data, (output / "multiqc_report.html").read_text(), output


@pytest.mark.parametrize("modules", [["kraken", "bracken"], ["bracken", "kraken"]])
def test_kraken_and_bracken_keep_distinct_values_for_the_same_sample(tmp_path, modules):
    inputs = tmp_path / "inputs"
    shutil.copytree(FIXTURES, inputs)
    # Original reports are excluded by GCF; the collector's .bracken copy is
    # authoritative here. A different sample name proves the ignored file is
    # not merely overwritten by the staged report.
    shutil.copy(inputs / "Synthetic.bracken", inputs / "Ignored.bracken_kreport")
    data, html, output = render(tmp_path, modules)
    raw = data["report_saved_raw_data"]
    assert set(raw["multiqc_kraken"]) == {"Synthetic"}
    assert set(raw["multiqc_bracken"]) == {"Synthetic"}
    assert raw["multiqc_kraken"]["Synthetic"]["S"]["Escherichia coli"] == 400
    assert raw["multiqc_bracken"]["Synthetic"]["S"]["Escherichia coli"] == 720
    stats = data["report_general_stats_data"]
    assert stats["kraken"]["Synthetic"]["pct_unclassified"] == 10
    assert stats["kraken"]["Synthetic"]["pct_top_one"] == 40
    assert stats["bracken"]["Synthetic"]["pct_top_one"] == 80
    assert data["report_general_stats_headers"]["bracken"]["pct_top_n"]["title"] == "Top 10 species"
    assert 'id="kraken-top-n"' not in html
    assert 'id="bracken-top-n"' in html
    assert (output / "multiqc_data" / "multiqc_kraken.txt").is_file()
    assert (output / "multiqc_data" / "multiqc_bracken.txt").is_file()
    columns = data["report_plot_data"]["general_stats_table"]["datasets"][0]["dt"]["section_by_id"]
    assert columns["bracken"]["column_by_key"]["pct_top_one"]["placement"] == 200
    assert columns["kraken"]["column_by_key"]["pct_top_one"]["hidden"] is True
    assert columns["kraken"]["column_by_key"]["pct_unclassified"]["placement"] == 220


def test_kraken_without_unclassified_is_not_misidentified_as_bracken(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    kraken = (FIXTURES / "Synthetic.kraken.kreport").read_text().splitlines(keepends=True)
    # All classified is valid Kraken output even though its first row is root.
    (inputs / "AllClassified.kraken.kreport").write_text("".join(kraken[1:]))
    shutil.copy(FIXTURES / "Synthetic.bracken", inputs)
    data, _, _ = render(tmp_path, ["kraken", "bracken"])
    raw = data["report_saved_raw_data"]
    assert set(raw["multiqc_kraken"]) == {"AllClassified"}
    assert set(raw["multiqc_bracken"]) == {"Synthetic"}


def test_bracken_selector_works_without_kraken(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    shutil.copy(FIXTURES / "Synthetic.bracken", inputs)
    data, html, _ = render(tmp_path, ["bracken"])
    assert "multiqc_bracken" in data["report_saved_raw_data"]
    assert "multiqc_kraken" not in data["report_saved_raw_data"]
    assert 'id="bracken"' in html


def test_search_pattern_override_preserves_distinct_modules(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    shutil.copy(FIXTURES / "Synthetic.kraken.kreport", inputs)
    shutil.copy(FIXTURES / "Synthetic.bracken", inputs / "Synthetic.estimated.kreport")
    data, _, _ = render(
        tmp_path,
        ["kraken", "bracken"],
        {
            "extra_fn_clean_exts": [".kraken.kreport", ".estimated.kreport"],
            # MultiQC 1.35 requires list form to override a list-valued default;
            # its config merge otherwise silently retains the default list.
            "sp": {"bracken": [{"fn": "*.estimated.kreport", "contents": "\troot", "num_lines": 1}]},
        },
    )
    raw = data["report_saved_raw_data"]
    assert set(raw["multiqc_kraken"]) == {"Synthetic"}
    assert set(raw["multiqc_bracken"]) == {"Synthetic"}
    assert raw["multiqc_kraken"]["Synthetic"]["S"]["Escherichia coli"] == 400
    assert raw["multiqc_bracken"]["Synthetic"]["S"]["Escherichia coli"] == 720
