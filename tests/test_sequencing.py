"""Run native early-report modules with GCF installed, without importing BFQ."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "sequencing"
RUN_ID = "260101_SYNTHETIC_0001_TEST"


def render(tmp_path, modules):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for module in modules:
        for fixture in (FIXTURES / module).iterdir():
            shutil.copy(fixture, inputs / fixture.name)
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
    return json.loads((out / "multiqc_data/multiqc_data.json").read_text())


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
