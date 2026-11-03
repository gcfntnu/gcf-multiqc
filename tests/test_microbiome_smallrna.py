"""Synthetic parser and installed-MultiQC regressions; no scientific workflows."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from gcf_multiqc.modules import qiime2, unitas

FIXTURES = Path(__file__).parent / "fixtures" / "microbiome_smallrna"


def test_dada2_region_identity_and_filtering_losses():
    data = qiime2.parse_dada2((FIXTURES / "qiime2/dada2/V3/metadata.tsv").read_text())
    assert set(data["V3"]) == {"001", "NA"}
    assert data["V3"]["001"] == {
        "filtered": 10,
        "denoised": 5,
        "merged": 5,
        "non-chimeric": 5,
        "passed": 75,
    }
    source = (FIXTURES / "qiime2/dada2/V4/metadata.tsv").read_text()
    with pytest.raises(ValueError, match="increase"):
        qiime2.parse_dada2(source.replace("50\t45", "40\t45"))
    with pytest.raises(ValueError, match="sample_region"):
        qiime2.parse_dada2(source.replace("001_V4", "001"))


def test_taxonomy_csv_quoting_metadata_and_duplicate_display_names():
    level, data, database = qiime2.parse_taxonomy((FIXTURES / "qiime2/taxa/level-2.csv").read_text(), "level-2.csv")
    assert level == "Phylum" and database == "greengenes"
    assert data["001"] == {"Alpha": 40, "Beta, subgroup": 30, "Unclassified": 5}
    assert "Sample_Group" not in data["NA"]
    _, genera, _ = qiime2.parse_taxonomy((FIXTURES / "qiime2/taxa/level-6.csv").read_text(), "level-6.csv")
    assert genera["001"] == {"Common": 75}


def test_ordination_coordinates_sample_identity_and_unbounded_groups():
    points = qiime2.parse_ordination((FIXTURES / "qiime2/deicode/ordination_samples.txt").read_text())
    assert points["001"]["x"] == -0.5
    assert points["001"]["name"] == "001"
    assert points["001"]["group"] == "group one"
    samples = [f"sample_{i:02}" for i in range(20)]
    coordinates = "\n".join(f"{sample}\t{i}\t{-i}" for i, sample in enumerate(samples))
    metadata = [f"{sample}\tgroup {i:02}" for i, sample in enumerate(samples)]
    prefix = f"Site\t20\t2\n{coordinates}\n\nSamples:\nSample_ID\tSample_Group\n"
    forward = qiime2.parse_ordination(prefix + "\n".join(metadata))
    reverse = qiime2.parse_ordination(prefix + "\n".join(reversed(metadata)))
    assert forward == reverse
    assert len({point["color"] for point in forward.values()}) == 20
    assert [point["name"] for point in forward.values()] == samples


@pytest.mark.parametrize(
    "contents, message",
    [
        ("Site\t1\t2\n001\t0\t1\n", "Samples:"),
        ("Site\t2\t2\n001\t0\t1\nSamples:\nSample_ID\n001\n", "Truncated"),
        ("Site\t1\t2\n001\tNaN\t1\nSamples:\nSample_ID\n001\n", "Non-finite"),
    ],
)
def test_malformed_ordination_fails_without_hanging(contents, message):
    with pytest.raises(ValueError, match=message):
        qiime2.parse_ordination(contents)


def test_unitas_per_sample_denominators_hierarchy_and_input_order():
    reports = {
        sample: (FIXTURES / f"unitas/{sample}/unitas.annotation_summary.txt").read_text() for sample in ("001", "1")
    }
    forward = {sample: unitas.parse_annotation(source) for sample, source in reports.items()}
    reverse = {sample: unitas.parse_annotation(source) for sample, source in reversed(list(reports.items()))}
    assert forward == reverse
    assert forward["001"] == ({"miRNA": 10, "rRNA": 90}, {"miRNA Percentage": 0.1})
    assert forward["1"] == ({"miRNA": 20, "rRNA": 150, "other": 30}, {"miRNA Percentage": 0.1})
    assert unitas.parse_annotation("miRNA\t0\nrRNA\t0\n")[1]["miRNA Percentage"] is None
    assert unitas.parse_annotation("rRNA\t10\n")[1]["miRNA Percentage"] == 0


def test_unitas_first_gene_threshold_and_sequence_counts():
    assert unitas.parse_mirna_counts((FIXTURES / "unitas/001/miR-table_synthetic.simplified.txt").read_text()) == 2
    assert unitas.parse_mirna_counts("miRNA\tCount\nsynthetic\t3\n") == 1
    assert unitas.parse_mirna_counts("miRNA\tCount\nsynthetic\t0\n") == 0
    assert unitas.parse_mirna_counts("miRNA\tCount\n") is None
    assert unitas.sequence_biotype("unitas.miR.synthetic.info") == "miRNA"
    assert unitas.sequence_biotype("unitas.non-miR.synthetic.info") is None
    assert unitas.parse_sequence_lengths(
        (FIXTURES / "unitas/001/sequence_info/unitas.miR.synthetic.info").read_text()
    ) == {21: 4, 22: 6}
    with pytest.raises(ValueError, match="nonnegative"):
        unitas.parse_mirna_counts("miRNA\tCount\nsynthetic\t-1\n")


def _run_multiqc(tmp_path, input_dir, modules, *extra_args):
    """Run the installed CLI with network/mail blocked and ambient config disabled."""
    guard = tmp_path / "guard"
    guard.mkdir()
    (guard / "sitecustomize.py").write_text(
        "import sys\n"
        "def protect(event, args):\n"
        "    if event in ('socket.connect', 'socket.getaddrinfo', 'smtplib.connect'):\n"
        "        raise RuntimeError('External network and SMTP are forbidden in synthetic tests')\n"
        "sys.addaudithook(protect)\n"
        "from multiqc import config\n"
        "config.find_user_files = lambda: None\n"
    )
    cache = tmp_path / "cache"
    cache.mkdir()
    environment = {key: value for key, value in os.environ.items() if not key.startswith("MULTIQC_")}
    environment.update(PYTHONPATH=str(guard), XDG_CACHE_HOME=str(cache), MPLCONFIGDIR=str(cache / "matplotlib"))
    command = [
        sys.executable,
        "-m",
        "multiqc",
        str(input_dir),
        "--strict",
        "--no-version-check",
        "--no-ai",
        "--no-ansi",
        "-k",
        "json",
        "-o",
        str(tmp_path / "out"),
    ]
    for module in modules:
        command.extend(["-m", module])
    command.extend(extra_args)
    result = subprocess.run(command, cwd=tmp_path, env=environment, text=True, capture_output=True, timeout=45)
    (tmp_path / "multiqc.log").write_text(result.stdout + result.stderr)
    return result, tmp_path / "out"


def _raw(out, filename):
    return json.loads((out / "multiqc_data" / "multiqc_data.json").read_text())["report_saved_raw_data"][filename]


def test_installed_report_contains_both_modules_and_correct_exports(tmp_path):
    result, out = _run_multiqc(tmp_path, FIXTURES, ["qiime2", "unitas"], "--require-logs")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Error validating config" not in result.stderr
    html = (out / "multiqc_report.html").read_text()
    for section in ("qiime2-dada2-stats", "qiime2-taxonomy", "qiime2-deicode", "unitas_seqlen", "unitas_annotation"):
        assert section in html
    dada2 = _raw(out, "multiqc_qiime2_dada2")
    assert dada2["Total"]["001"]["passed"] == 105
    assert dada2["Total"]["NA"]["passed"] == 120
    assert _raw(out, "multiqc_qiime2_taxonomy")["Phylum"]["001"]["Beta, subgroup"] == 30
    assert _raw(out, "multiqc_qiime2_rpca")["001"]["name"] == "001"
    metrics = _raw(out, "multiqc_unitas")
    assert set(metrics) == {"001", "1"}
    assert metrics["001"] == {"miRNA Percentage": 0.1, "miRNA Detected": 2}
    assert metrics["1"] == {"miRNA Percentage": 0.1, "miRNA Detected": 1}
    report = json.loads((out / "multiqc_data/multiqc_data.json").read_text())
    assert report["report_plot_data"]["dada2-rpca-plot"]["layout"]["showlegend"] is True
    points = report["report_plot_data"]["dada2-rpca-plot"]["datasets"][0]["points"]
    assert {point["name"] for point in points} == {"001", "NA"}
    assert {point["group"] for point in points} == {"group one", "group two"}
    assert all(point["color"] == qiime2.group_color(point["group"]) for point in points)


def test_report_keeps_available_sections_and_zero_is_unavailable(tmp_path):
    inputs = tmp_path / "inputs"
    sample = inputs / "unitas" / "NA"
    sample.mkdir(parents=True)
    (sample / "unitas.annotation_summary.txt").write_text("miRNA\t0\nrRNA\t0\n")
    taxonomy = inputs / "taxa"
    taxonomy.mkdir()
    shutil.copy(FIXTURES / "qiime2/taxa/level-2.csv", taxonomy / "level-2.csv")
    result, out = _run_multiqc(tmp_path, inputs, ["qiime2", "unitas"], "--require-logs")
    assert result.returncode == 0, result.stdout + result.stderr
    assert _raw(out, "multiqc_unitas")["NA"]["miRNA Percentage"] is None
    assert _raw(out, "multiqc_qiime2_taxonomy")["Phylum"]["001"]["Alpha"] == 40
    raw = json.loads((out / "multiqc_data/multiqc_data.json").read_text())["report_saved_raw_data"]
    assert "multiqc_qiime2_rpca" not in raw


def test_report_sample_ignore_applies_to_every_nested_dataset(tmp_path):
    result, out = _run_multiqc(tmp_path, FIXTURES, ["qiime2", "unitas"], "--ignore-samples", "001")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "001" not in _raw(out, "multiqc_unitas")
    assert "001" not in _raw(out, "multiqc_qiime2_rpca")
    for samples in _raw(out, "multiqc_unitas_lengths").values():
        assert "001" not in samples
    for samples in _raw(out, "multiqc_qiime2_dada2").values():
        assert "001" not in samples
    for samples in _raw(out, "multiqc_qiime2_taxonomy").values():
        assert "001" not in samples


def test_malformed_ordination_installed_cli_fails_promptly(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "ordination_samples.txt").write_text("Site\t1\t2\n001\t0\t1\n")
    result, out = _run_multiqc(tmp_path, inputs, ["qiime2"])
    assert result.returncode != 0
    assert "missing the appended Samples:" in result.stderr
    assert not (out / "multiqc_report.html").exists()


@pytest.mark.parametrize("module", ["qiime2", "unitas"])
def test_empty_inputs_are_cleanly_no_samples(tmp_path, module):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    result, out = _run_multiqc(tmp_path, inputs, [module], "--require-logs")
    assert result.returncode != 0
    assert "no log files were found" in result.stderr.lower()
    assert "Traceback" not in result.stderr
    assert not (out / "multiqc_report.html").exists()
