"""Render public synthetic inputs through installed MultiQC and GCF entry points."""

from copy import deepcopy
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
        *[str(path) for path in (inputs if isinstance(inputs, list) else [inputs])],
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
    assert data["report_saved_raw_data"]["multiqc_fastqc"]["Synthetic_R1"]["Total Sequences"] == 1000
    assert (output / "multiqc_data" / "multiqc_fastqc.txt").is_file()
    assert data["report_plot_data"]


def paired_fastqc_archives(tmp_path):
    archives = []
    for mate, mean_quality, status in ((1, 35, "warn"), (2, 25, "pass")):
        name = f"Synthetic_R{mate}"
        content = (FIXTURES / "fastqc/fastqc_data.txt").read_text()
        content = content.replace("Synthetic_R1", name)
        if mate == 2:
            content = content.replace("1000", "2000").replace("75000", "150000")
        content += (
            ">>Per base sequence quality\tpass\n"
            "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n"
            f"1\t{mean_quality}\t{mean_quality}\t20\t40\t15\t40\n"
            ">>END_MODULE\n"
            f">>Overrepresented sequences\t{status}\n"
        )
        if mate == 1:
            content += "#Sequence\tCount\tPercentage\tPossible Source\nACGTACGTACGT\t200\t20.0\tSynthetic control\n"
        content += ">>END_MODULE\n"
        archive_path = tmp_path / f"{name}_fastqc.zip"
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr(f"{name}_fastqc/", "")
            archive.writestr(f"{name}_fastqc/fastqc_data.txt", content)
        archives.append(archive_path)
    return archives


@pytest.mark.parametrize("reverse", [False, True])
def test_fastqc_alias_keeps_both_mates_independent_of_input_order(tmp_path, reverse):
    archives = paired_fastqc_archives(tmp_path)
    originals = {path: path.read_bytes() for path in archives}
    inputs = list(reversed(archives)) if reverse else archives
    data, html, _ = render(
        tmp_path,
        [*inputs, FIXTURES / "Synthetic.fastp.json"],
        ["fastqc_rnaseq", "fastp"],
        {"module_order": ["fastqc_rnaseq", "fastp"]},
    )
    raw = data["report_saved_raw_data"]["multiqc_fastqc"]
    assert {sample: values["Total Sequences"] for sample, values in raw.items()} == {
        "Synthetic_R1": 1000,
        "Synthetic_R2": 2000,
    }
    assert raw["Synthetic_R1"]["overrepresented_sequences"] == "warn"
    assert raw["Synthetic_R2"]["overrepresented_sequences"] == "pass"
    assert set(data["report_general_stats_data"]) == {"fastp"}
    assert set(data["report_general_stats_data"]["fastp"]) == {"Synthetic"}
    quality = data["report_plot_data"]["fastqc_rnaseq_per_base_sequence_quality_plot"]["datasets"][0]
    assert {line["name"]: line["pairs"] for line in quality["lines"]} == {
        "Synthetic_R1": [[1, 35.0]],
        "Synthetic_R2": [[1, 25.0]],
    }
    overrep = data["report_plot_data"]["fastqc_rnaseq_overrepresented_sequences_plot"]["datasets"][0]
    top = next(category for category in overrep["cats"] if category["name"] == "Top overrepresented sequence")
    assert dict(zip(overrep["samples"], top["data"])) == {"Synthetic_R1": 20, "Synthetic_R2": 0}
    assert 'id="fastqc_rnaseq_top_overrepresented_sequences"' in html
    assert {path: path.read_bytes() for path in archives} == originals


def test_native_fastqc_optout_keeps_native_mate_names_and_general_stats(tmp_path):
    data, html, _ = render(tmp_path, paired_fastqc_archives(tmp_path), ["fastqc"], {"gcf_compatibility": False})
    assert set(data["report_saved_raw_data"]["multiqc_fastqc"]) == {"Synthetic_R1", "Synthetic_R2"}
    assert set(data["report_general_stats_data"]["fastqc"]) == {"Synthetic_R1", "Synthetic_R2"}
    assert 'id="fastqc"' in html
    assert 'id="fastqc_rnaseq"' not in html


def test_fastqc_alias_does_not_mutate_global_name_cleaning(monkeypatch):
    from multiqc import config, report
    from multiqc.base_module import BaseMultiqcModule

    from gcf_multiqc.modules.fastqc_rnaseq import MultiqcModule

    rules = [".gz", ".fastq", "_R1", "_R2", {"type": "truncate", "pattern": ".trimmed", "module": "fastqc_rnaseq"}]
    original = deepcopy(rules)
    monkeypatch.setattr(config, "fn_clean_exts", rules)
    monkeypatch.setattr(config, "fn_clean_sample_names", True)
    monkeypatch.setattr(config, "fn_clean_trim", [])
    monkeypatch.setattr(config, "prepend_dirs", False)
    monkeypatch.setattr(config, "use_filename_as_sample_name", False)
    monkeypatch.setattr(config, "sample_names_replace", {})
    monkeypatch.setattr(report, "sample_names", [])
    module = MultiqcModule.__new__(MultiqcModule)
    module.anchor = "fastqc_rnaseq"
    source = {"root": "/synthetic", "fn": "Synthetic_R1_fastqc.zip", "sp_key": "fastqc_rnaseq/zip"}
    assert module.clean_s_name("Synthetic_R1.fastq.gz", source) == "Synthetic_R1"
    assert module.clean_s_name("Synthetic_R2.fastq.gz", source) == "Synthetic_R2"
    assert config.fn_clean_exts is rules and config.fn_clean_exts == original
    standard = BaseMultiqcModule.__new__(BaseMultiqcModule)
    standard.anchor = "fastp"
    assert standard.clean_s_name("Synthetic_R1.fastq.gz", source) == "Synthetic"


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
