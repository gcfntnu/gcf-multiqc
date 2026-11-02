"""Run native early-report modules with GCF installed, without importing BFQ."""

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "sequencing"
RUN_ID = "260101_SYNTHETIC_0001_TEST"


def render(tmp_path, modules, prepare=None):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for module in modules:
        for fixture in (FIXTURES / module).iterdir():
            shutil.copy(fixture, inputs / fixture.name)
    if prepare:
        prepare(inputs)
    original_inputs = {path: path.read_bytes() for path in inputs.rglob("*") if path.is_file()}
    # The inherited test guard also protects children. Disable ambient personal
    # config to keep report fixtures independent of the developer's account.
    bootstrap = (
        "import sys, runpy\n"
        "def guard(event, args):\n"
        "    if event in {'socket.connect', 'socket.getaddrinfo', 'smtplib.connect'}:\n"
        "        raise RuntimeError('Network and SMTP disabled in synthetic report tests')\n"
        "sys.addaudithook(guard)\n"
        "from multiqc import config\n"
        "config.find_user_files = lambda: None\n"
        "runpy.run_module('multiqc', run_name='__main__')\n"
    )
    out = tmp_path / "out"
    command = [
        sys.executable,
        "-c",
        bootstrap,
        "--strict",
        "--interactive",
        "--no-version-check",
        "--no-megaqc-upload",
        "--no-ai",
        "--no-ansi",
        "--require-logs",
        "-o",
        str(out),
    ]
    for module in modules:
        command += ["-m", module]
    command.append(str(inputs))
    environment = {key: value for key, value in os.environ.items() if not key.startswith("MULTIQC_")}
    environment.update(XDG_CACHE_HOME=str(tmp_path / "cache"), MPLCONFIGDIR=str(tmp_path / "mpl"))
    result = subprocess.run(command, cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=45)
    (tmp_path / "multiqc.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "GCF compatibility enabled" in result.stderr
    assert "Error validating config" not in result.stderr
    assert (out / "multiqc_report.html").stat().st_size > 1000
    assert {path: path.read_bytes() for path in inputs.rglob("*") if path.is_file()} == original_inputs
    return json.loads((out / "multiqc_data/multiqc_data.json").read_text())


def write_indexed_demux(inputs, rows, quality="modern", run_id=RUN_ID):
    inputs.mkdir(exist_ok=True)
    (inputs / "RunInfo.xml").write_text((FIXTURES / "bclconvert/RunInfo.xml").read_text().replace(RUN_ID, run_id))
    demux_rows, quality_rows = [], []
    for lane, sample, index, count, perfect, mean_q in rows:
        demux = {
            "Lane": lane,
            "SampleID": sample,
            "Index": index,
            "# Reads": count,
            "# Perfect Index Reads": perfect,
            "# One Mismatch Index Reads": count - perfect,
        }
        if quality == "embedded":
            demux.update({"# of >= Q30 Bases (PF)": count * 240, "Mean Quality Score (PF)": mean_q})
        demux_rows.append(demux)
        quality_rows.append(
            {
                "Lane": lane,
                "SampleID": sample,
                "Index": index,
                "Yield": count * 302,
                "YieldQ30": count * 240,
                "QualityScoreSum": count * 302 * mean_q,
            }
        )
    files = {"Demultiplex_Stats.csv": demux_rows}
    if quality == "modern":
        files["Quality_Metrics.csv"] = quality_rows
    for name, values in files.items():
        with (inputs / name).open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(values[0]))
            writer.writeheader()
            writer.writerows(values)


@pytest.mark.parametrize("quality", ["modern", "embedded", "missing"])
def test_bclconvert_sums_repeated_indexes_in_lane_plots_and_exports(tmp_path, quality):
    indexes = ["AAAA+CCCC", "AAAC+CCCC", "AAAG+CCCC", "AAAT+CCCC"]
    rows = [(1, "pooled", index, 10 * i, 8 * i, 10 * (i + 1)) for i, index in enumerate(indexes, 1)]
    rows += [
        (2, "pooled", indexes[0], 5, 4, 20),
        (2, "pooled", indexes[1], 15, 12, 40),
        (1, "single", "GGGG+TTTT", 8, 8, 30),
        (1, "Undetermined", "", 12, 0, 0),
        (2, "Undetermined", "", 30, 0, 0),
    ]
    report = render(tmp_path, ["bclconvert"], lambda inputs: write_indexed_demux(inputs, rows, quality))
    raw = report["report_saved_raw_data"]
    pooled = raw["multiqc_bclconvert_bysample"]["pooled"]
    assert pooled["clusters"] == 120  # clusters, not twice that number of paired-end reads
    assert pooled["yield_"] == 120 * 302
    assert set(pooled["index"].split(", ")) == set(indexes)
    for lane, count, mean_q in (("L1", 100, 40), ("L2", 20, 35)):
        chunk = pooled["lanes"][lane]
        assert chunk == raw["multiqc_bclconvert_bylane"][f"{RUN_ID} - {lane}"]["samples"]["pooled"]
        assert chunk["clusters"] == count
        assert chunk["yield_"] == count * 302
        assert chunk["perfect_index_reads"] == count * 0.8
        assert chunk["one_mismatch_index_reads"] == count * 0.2
        assert chunk["percent_perfect_index_reads"] == 80
        assert chunk["percent_one_mismatch_index_reads"] == 20
        if quality == "missing":
            assert chunk["yield_q30"] is None and chunk["mean_quality"] is None
        else:
            assert chunk["yield_q30"] == count * 240
            assert chunk["percent_yield_q30"] == pytest.approx(240 / 302 * 100)
            assert chunk["mean_quality"] == mean_q
    if quality != "missing":
        assert pooled["mean_quality"] == pytest.approx((100 * 40 + 20 * 35) / 120)
    per_lane = report["report_plot_data"]["bclconvert_sample_counts"]["datasets"][1]
    plotted = {category["name"]: dict(zip(per_lane["samples"], category["data"])) for category in per_lane["cats"]}
    assert plotted["L1"]["pooled"] == 100
    assert plotted["L2"]["pooled"] == 20
    assert plotted["L1"]["single"] == 8
    lane_counts = report["report_plot_data"]["bclconvert_lane_counts"]["datasets"][0]
    undetermined = next(category for category in lane_counts["cats"] if category["name"] == "Undetermined Reads")
    assert dict(zip(lane_counts["samples"], undetermined["data"])) == {"L1": 12, "L2": 30}


@pytest.mark.parametrize("same_run", [True, False])
def test_bclconvert_keeps_chunks_distinct_across_files_and_runs(tmp_path, same_run):
    second_run = RUN_ID if same_run else "260102_SYNTHETIC_0002_TEST"

    def prepare(inputs):
        write_indexed_demux(inputs, [(1, "pooled", "AAAA+CCCC", 10, 8, 20), (1, "Undetermined", "", 90, 0, 0)])
        write_indexed_demux(
            inputs / "second",
            [(1, "pooled", "GGGG+TTTT", 30, 24, 40), (1, "Undetermined", "", 70, 0, 0)],
            run_id=second_run,
        )

    report = render(tmp_path, ["bclconvert"], prepare)
    raw = report["report_saved_raw_data"]
    pooled = raw["multiqc_bclconvert_bysample"]["pooled"]
    assert pooled["clusters"] == 40 and pooled["yield_"] == 40 * 302
    assert pooled["mean_quality"] == 35
    assert set(pooled["index"].split(", ")) == {"AAAA+CCCC", "GGGG+TTTT"}
    expected = {"L1": 40} if same_run else {f"{RUN_ID} - L1": 10, f"{second_run} - L1": 30}
    assert {lane: chunk["clusters"] for lane, chunk in pooled["lanes"].items()} == expected
    for lane, count in expected.items():
        assert pooled["lanes"][lane]["yield_"] == count * 302
        assert pooled["lanes"][lane]["yield_q30"] == count * 240
    per_lane = report["report_plot_data"]["bclconvert_sample_counts"]["datasets"][1]
    assert per_lane["samples"] == ["pooled"]
    assert {category["name"]: category["data"][0] for category in per_lane["cats"]} == expected


def test_native_bclconvert_preserves_zero_and_missing_quality(tmp_path):
    report = render(tmp_path, ["bclconvert"])
    raw = report["report_saved_raw_data"]
    samples = raw["multiqc_bclconvert_bysample"]
    assert set(samples) == {"s1", "s2"}
    assert samples["s1"]["clusters"] == 1
    assert samples["s1"]["perfect_index_reads"] == 1
    assert samples["s1"]["cluster_length"] == 302
    assert samples["s1"]["yield_"] == 302
    assert samples["s2"]["clusters"] == 0
    assert samples["s2"]["yield_"] == 0
    assert samples["s1"]["yield_q30"] is None
    assert samples["s1"]["percent_yield_q30"] is None
    lane = raw["multiqc_bclconvert_bylane"][f"{RUN_ID} - L1"]
    assert lane["clusters"] == 1
    assert lane["top_unknown_barcodes"] == {"AAAA-GGGG": 9000}
    categories = {
        category["name"]: category
        for category in report["report_plot_data"]["bclconvert_lane_counts"]["datasets"][0]["cats"]
    }
    assert categories["Undetermined Reads"]["data"] == [9999]
    assert categories["Undetermined Reads"]["data_pct"] == [99.99]


def test_native_bcl2fastq_retains_assigned_and_undetermined_counts(tmp_path):
    report = render(tmp_path, ["bcl2fastq"])
    raw = report["report_saved_raw_data"]
    samples = raw["multiqc_bcl2fastq_bysample"]
    assert samples["s1-first"]["total"] == 75
    assert samples["s1-first"]["mean_qscore"] == 35
    assert samples["s1-first"]["percent_Q30"] == pytest.approx(140 / 151 * 100)
    assert samples["s2-second"]["total"] == 0
    assert samples["s2-second"]["mean_qscore"] == "NA"
    assert samples["undetermined"]["total"] == 25
    lane = raw["multiqc_bcl2fastq_bylane"][f"{RUN_ID} - L1"]
    assert lane["total"] == 75
    assert lane["undetermined"] == 25
    assert lane["unknown_barcodes"] == {"ATGC+ATGC": 20}


def test_native_interop_uses_aggregate_lane_metrics(tmp_path):
    report = render(tmp_path, ["interop"])
    summary = report["report_saved_raw_data"]["interop_runsummary"]["interop"]
    assert summary["summary"]["Non-indexed"]["Yield"] == 30.2
    assert summary["summary"]["Non-indexed"]["%>=Q30"] == 93.5
    assert summary["summary"]["Total"]["Yield"] == 31.8
    assert len(summary["details"]) == 4
    lane = summary["details"]["Lane 1 - Read 1"]
    assert lane["Reads"] == 100
    assert lane["%>=Q30"] == 94  # The surface-specific fixture row has Q30=20.


@pytest.mark.parametrize("demultiplexer", ["bclconvert", "bcl2fastq"])
def test_native_early_modules_render_together(tmp_path, demultiplexer):
    report = render(tmp_path, [demultiplexer, "interop"])
    assert "interop_runsummary" in report["report_saved_raw_data"]
    assert f"multiqc_{demultiplexer}_bysample" in report["report_saved_raw_data"]
