"""Large reports use flat plots; an interactive-only smoke misses exporter failures."""

from pathlib import Path
import subprocess
import sys


def test_flat_report_embeds_real_images(tmp_path):
    fixtures = Path(__file__).parent / "fixtures" / "compatibility"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "multiqc",
            str(fixtures),
            "-m",
            "fastp",
            "--strict",
            "--require-logs",
            "--flat",
            "--no-version-check",
            "--no-ai",
            "--no-megaqc-upload",
            "-o",
            str(tmp_path / "report"),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    html = (tmp_path / "report/multiqc_report.html").read_text()
    assert "data:image/png;base64," in html
    assert "Unable to export" not in result.stderr
