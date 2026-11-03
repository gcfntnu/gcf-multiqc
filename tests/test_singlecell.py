"""Synthetic single-cell adapter contracts, including real installed CLI rendering."""

import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import Mock

import pytest

from gcf_multiqc.modules import cellranger_count, cellranger_mkfastq, parse, starsolo

FIXTURES = Path(__file__).parent / "fixtures" / "singlecell"


@pytest.mark.parametrize("filename", ["parse/001/agg_sample_summary.csv", "parse/NA.agg_samp_ana_summary.csv"])
def test_parse_both_csv_orientations(filename):
    with (FIXTURES / filename).open() as handle:
        data = parse.parse_summary(handle)
    assert data["number_of_reads"] == 1000
    assert data["number_of_cells"] == 10
    assert data["sequencing_saturation"] == 0.25
    assert data["transcriptome_map_fraction"] == 0.85


@pytest.mark.parametrize(
    "filename,root,expected",
    [
        ("agg_sample_summary.csv", "/synthetic/001", "001"),
        ("agg_samp_ana_summary.csv", "/synthetic/NA", "NA"),
        ("001.agg_samp_ana_summary.csv", "/synthetic/wrong", "001"),
        ("NA.agg_sample_summary.csv", "/synthetic/wrong", "NA"),
        ("all_samples.agg_sample_summary.csv", "/synthetic/logs", "all_samples"),
    ],
)
def test_parse_identity_is_literal_and_filename_prefix_wins(filename, root, expected):
    assert parse.sample_name(filename, root) == expected


@pytest.mark.parametrize(
    "contents,match",
    [
        ("sample,number_of_reads\nS1,10\n", "all-sample"),
        ("sample,number_of_reads\nall-sample,1\nall-sample,2\n", "duplicate"),
        ("sample,number_of_reads,number_of_reads\nall-sample,1,2\n", "duplicate"),
        ("metric,all-sample\nnumber_of_reads,1\nnumber_of_reads,2\n", "duplicate"),
        ("sample,number_of_reads\nall-sample,1,2\n", "column counts"),
        ("sample,number_of_reads\nall-sample,unknown\n", "Invalid Parse metric"),
        ("sample,number_of_reads,sequencing_saturation\nall-sample,1,25\n", "Out-of-range"),
    ],
)
def test_parse_rejects_ambiguous_or_malformed_csv(contents, match):
    with pytest.raises(ValueError, match=match):
        parse.parse_summary(io.StringIO(contents))


def test_parse_fraction_and_missing_value_semantics():
    data = parse.parse_summary(
        io.StringIO(
            "sample,number_of_reads,number_of_cells,sequencing_saturation,bc1_Q30,annotation\n"
            'all-sample,"1,000",0,25%,N/A,<1\n'
        )
    )
    assert data == {"number_of_reads": 1000, "number_of_cells": 0, "sequencing_saturation": 0.25, "annotation": "<1"}


def test_parse_splitpipe_producer_columns_preserve_extra_metrics_and_all_sample_summary():
    # Column names follow the collected SplitPipe format; every value is invented.
    metrics = {
        "sample_well_count": 2,
        "number_of_cells": 10,
        "GRCh38_median_tscp_per_cell": 200,
        "GRCh38_median_tscp_at50": 250,
        "GRCh38_median_genes_per_cell": 100,
        "mean_reads_per_cell": 100,
        "number_of_reads": 1000,
        "number_of_tscp": 500,
        "sequencing_saturation": 0.25,
        "valid_barcode_fraction": 0.95,
        "bc1_Q30": 0.96,
        "bc2_Q30": 0.97,
        "bc3_Q30": 0.98,
        "cDNA_Q30": 0.99,
        "polyN_Q30": 0.90,
        "sample_tmap_fraction": 0.75,
        "transcriptome_map_fraction": 0.85,
        "tso_fraction_in_read1": 0.05,
        "GRCh38_fraction_reads_in_cells": 0.80,
        "GRCh38_fraction_tscp_in_cells": 0.82,
        "GRCh38_fraction_exonic": 0.72,
        "cell_tscp_cutoff": 50,
        "cell_tscp_f01_slope": -0.75,
    }
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(["sample", *metrics])
    # The aggregate can follow biological-sample rows; it must be the selected row.
    writer.writerow(["001", *([0] * len(metrics))])
    writer.writerow(["all-sample", *metrics.values()])
    stream.seek(0)
    assert parse.parse_summary(stream) == metrics


def test_starsolo_numeric_metrics_and_legacy_spelling():
    with (FIXTURES / "starsolo/001_Summary.csv").open() as handle:
        data = starsolo.parse_summary(handle)
    assert data["Mean Reads per Cell"] == 85.5
    assert data["Sequencing Saturation"] == 0.25
    assert data["Reads Mapped to Transcriptome: Unique+Multipe Genes"] == 0.85
    assert starsolo.sample_name("001_Summary.csv", "/wrong") == "001"
    assert starsolo.sample_name("Summary.csv", "/synthetic/NA") == "NA"


@pytest.mark.parametrize(
    "contents,match",
    [
        ("Number of Reads,10\nNumber of Reads,20\n", "Duplicate"),
        ("Number of Reads,10,20\n", "metric and value"),
        ("Number of Reads,NaN\n", "missing Number"),
        ("Number of Reads,inf\n", "Out-of-range"),
        ("Number of Reads,10\nSequencing Saturation,25\n", "Out-of-range"),
    ],
)
def test_starsolo_rejects_invalid_metrics(contents, match):
    with pytest.raises(ValueError, match=match):
        starsolo.parse_summary(io.StringIO(contents))


def test_count_normalizes_numbers_and_percentages_without_guessing_fractions():
    with (FIXTURES / "count/001.metrics_summary.csv").open() as handle:
        data = cellranger_count.parse_summary(handle)
    assert data["Number of Reads"] == 1000
    assert data["Median UMI Counts per Cell"] == 1234
    assert data["Valid Barcodes"] == 98
    assert data["Sequencing Saturation"] == 25
    assert data["Reads Mapped Confidently to Intergenic Regions"] == 5
    # A bare 0.5 in a percent column is 0.5%, not an inferred fraction.
    data = cellranger_count.parse_summary(io.StringIO("Estimated Number of Cells,Valid Barcodes\n0,0.5\n"))
    assert data == {"Estimated Number of Cells": 0, "Valid Barcodes": 0.5}


@pytest.mark.parametrize(
    "filename,root,expected",
    [
        ("001.metrics_summary.csv", "/synthetic/wrong", "001"),
        ("NA.metrics_summary.csv", "/synthetic/wrong", "NA"),
        ("metrics_summary.csv", "/synthetic/NULL/outs", "NULL"),
        ("metrics_summary.csv", "/synthetic/001", "001"),
    ],
)
def test_count_sample_identity(filename, root, expected):
    assert cellranger_count.sample_name(filename, root) == expected


@pytest.mark.parametrize(
    "contents,match",
    [
        ("Estimated Number of Cells\n1\n2\n", "exactly one"),
        ("Estimated Number of Cells,Estimated Number of Cells\n1,2\n", "duplicate"),
        ("Estimated Number of Cells\n1,2\n", "column counts"),
        ("Estimated Number of Cells,Valid Barcodes\n1,101%\n", "Out-of-range"),
        ("Estimated Number of Cells\nnot measured\n", "Invalid Cell Ranger"),
    ],
)
def test_count_rejects_malformed_or_multiple_records(contents, match):
    with pytest.raises(ValueError, match=match):
        cellranger_count.parse_summary(io.StringIO(contents))


def test_count_missing_optional_regions_are_not_invented_as_zero(caplog):
    module = cellranger_count.MultiqcModule.__new__(cellranger_count.MultiqcModule)
    module.cellranger_count_data = {"001": {"Estimated Number of Cells": 0}}
    assert module.cellranger_count_mapped_chart() is None
    assert "missing region metrics: 001" in caplog.text
    assert module.cellranger_count_data == {"001": {"Estimated Number of Cells": 0}}


def test_mkfastq_preserves_json_identifiers_fractions_and_zero():
    with (FIXTURES / "mkfastq/qc_summary.json").open() as handle:
        version, data = cellranger_mkfastq.parse_summary(handle)
    assert version == "cellranger-atac-2.0.0"
    assert set(data) == {"001", "NA"}
    assert data["001"]["barcode_exact_match_ratio"] == 0.98
    assert data["NA"]["number_reads"] == data["NA"]["gem_count_estimate"] == 0


@pytest.mark.parametrize(
    "payload,match",
    [
        ({"sample_qc": {}}, "10x_software_version"),
        ({"10x_software_version": "fixture", "sample_qc": {"S1": {}}}, "all metrics"),
        ({"10x_software_version": "fixture", "sample_qc": {"S1": {"all": {}}}}, "number_reads"),
        ({"10x_software_version": "fixture", "sample_qc": {"S1": {"all": {"number_reads": True}}}}, "Invalid mkfastq"),
        ({"10x_software_version": "fixture", "sample_qc": {"S1": {"all": {"number_reads": -1}}}}, "Out-of-range"),
    ],
)
def test_mkfastq_rejects_incomplete_or_invalid_objects(payload, match):
    with pytest.raises(ValueError, match=match):
        cellranger_mkfastq.parse_summary(io.StringIO(json.dumps(payload)))


def test_mkfastq_rejects_duplicate_json_identifiers_before_overwrite():
    with pytest.raises(ValueError, match="Duplicate mkfastq JSON key '001'"):
        cellranger_mkfastq.parse_summary(
            io.StringIO(
                '{"10x_software_version":"fixture","sample_qc":{"001":{"all":{"number_reads":1}},'
                '"001":{"all":{"number_reads":2}}}}'
            )
        )


@pytest.mark.parametrize(
    "adapter,attribute,method,fixture,filename",
    [
        (parse, "parse_data", "parse_log", "parse/001/agg_sample_summary.csv", "001.agg_sample_summary.csv"),
        (starsolo, "starsolo_data", "parse_starsolo_log", "starsolo/001_Summary.csv", "001_Summary.csv"),
        (
            cellranger_count,
            "cellranger_count_data",
            "parse_cellranger_count_log",
            "count/001.metrics_summary.csv",
            "001.metrics_summary.csv",
        ),
        (cellranger_mkfastq, "cellranger_data", "parse_cellranger_log", "mkfastq/qc_summary.json", "qc_summary.json"),
    ],
)
def test_duplicate_inputs_never_overwrite(adapter, attribute, method, fixture, filename):
    module = adapter.MultiqcModule.__new__(adapter.MultiqcModule)
    setattr(module, attribute, {})
    module._sources = {}
    module.cellranger_versions = set()
    module.add_data_source = Mock()
    module.add_software_version = Mock()
    content = (FIXTURES / fixture).read_text()
    first = {"fn": filename, "root": "/synthetic/first", "f": io.StringIO(content)}
    getattr(module, method)(first)
    original = json.loads(json.dumps(getattr(module, attribute)))
    with pytest.raises(ValueError, match="Duplicate.*001") as error:
        getattr(module, method)({"fn": filename, "root": "/synthetic/second", "f": io.StringIO(content)})
    assert "first" in str(error.value) and "second" in str(error.value)
    assert getattr(module, attribute) == original


def _run_multiqc(tmp_path, inputs, modules, extra_config=""):
    # Child processes cannot reach SMTP, version servers or any other network.
    guard = tmp_path / "guard"
    guard.mkdir()
    (guard / "sitecustomize.py").write_text(
        "import sys\n"
        "def guard(event, args):\n"
        "    if event in ('socket.connect', 'smtplib.connect'):\n"
        "        raise RuntimeError('Network disabled in synthetic MultiQC tests')\n"
        "sys.addaudithook(guard)\n"
    )
    config = tmp_path / "config.yaml"
    config.write_text(
        "no_version_check: true\nmegaqc_upload: false\nread_count_multiplier: 1\nread_count_prefix: ''\n" + extra_config
    )
    env = os.environ.copy()
    env.update({"HOME": str(tmp_path), "XDG_CACHE_HOME": str(tmp_path / "cache"), "PYTHONPATH": str(guard)})
    command = [
        str(Path(sys.executable).with_name("multiqc")),
        str(inputs),
        "--config",
        str(config),
        "--outdir",
        str(tmp_path / "out"),
        "--strict",
        "--require-logs",
        "--interactive",
        "--no-megaqc-upload",
    ]
    for module in modules:
        command.extend(["-m", module])
    result = subprocess.run(command, cwd=tmp_path, env=env, text=True, capture_output=True, timeout=90)
    (tmp_path / "multiqc.log").write_text(result.stdout + result.stderr)
    return result


@pytest.mark.parametrize("mkfastq_selector", ["cellranger", "gcf_cellranger_mkfastq"])
def test_installed_cli_renders_legacy_singlecell_contracts(tmp_path, mkfastq_selector):
    result = _run_multiqc(tmp_path, FIXTURES, ["parse", "starsolo", "cellranger_count", mkfastq_selector])
    assert result.returncode == 0, result.stdout + result.stderr
    html = (tmp_path / "out/multiqc_report.html").read_text()
    for anchor in (
        "parse-summary",
        "starsolo-summary",
        "cellranger_count-summary",
        "cellranger_count_confidently_mapped",
        f"{mkfastq_selector}-mkfastq-qc",
    ):
        assert anchor in html
    data = json.loads((tmp_path / "out/multiqc_data/multiqc_data.json").read_text())
    raw = data["report_saved_raw_data"]
    assert set(raw["multiqc_parse"]) == {"001", "NA", "all_samples"}
    assert raw["multiqc_parse"]["all_samples"]["number_of_reads"] == 2000
    assert raw["multiqc_starsolo"]["001"]["Mean Reads per Cell"] == 85.5
    assert raw["multiqc_cellranger_count"]["001"]["Number of Reads"] == 1000
    assert raw["multiqc_cellranger_count"]["001"]["Valid Barcodes"] == 98
    assert raw[f"multiqc_{mkfastq_selector}"]["NA"]["gem_count_estimate"] == 0
    assert raw[f"multiqc_{mkfastq_selector}"]["001"]["barcode_exact_match_ratio"] == 0.98
    assert "GEM estimate" in html  # Metric absent from first sample must still be shown.


def test_explicit_mkfastq_alias_coexists_with_upstream_cellranger(tmp_path):
    result = _run_multiqc(
        tmp_path, FIXTURES, ["cellranger", "gcf_cellranger_mkfastq"], extra_config="gcf_compatibility: false\n"
    )
    assert result.returncode == 0, result.stdout + result.stderr
    html = (tmp_path / "out/multiqc_report.html").read_text()
    assert "cellranger-count-stats" in html
    assert "gcf_cellranger_mkfastq-mkfastq-qc" in html
    data = json.loads((tmp_path / "out/multiqc_data/multiqc_data.json").read_text())["report_saved_raw_data"]
    assert data["multiqc_gcf_cellranger_mkfastq"]["001"]["number_reads"] == 1000
    assert data["multiqc_cellranger_count"]["SYNTHETIC_NATIVE"]["reads"] == 2000


def test_installed_cli_fails_on_malformed_parse_instead_of_silent_omission(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "S1.agg_sample_summary.csv").write_text("sample,number_of_reads\nmissing-all-sample,10\n")
    result = _run_multiqc(tmp_path, inputs, ["parse"])
    assert result.returncode != 0
    assert "S1.agg_sample_summary.csv" in result.stderr
    assert "all-sample" in result.stderr
