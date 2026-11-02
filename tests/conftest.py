"""Keep all local report tests offline, including their Python subprocesses."""
import os
from pathlib import Path
import sys

import pytest

SUPPORT = Path(__file__).parent / "support"
sys.path.insert(0, str(SUPPORT))
import sitecustomize  # noqa: E402, F401


@pytest.fixture(autouse=True)
def offline_children(monkeypatch, tmp_path):
    inherited = os.environ.get("PYTHONPATH", "")
    monkeypatch.setenv("PYTHONPATH", str(SUPPORT) + (os.pathsep + inherited if inherited else ""))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mpl"))
    monkeypatch.setenv("MULTIQC_NO_VERSION_CHECK", "true")
    monkeypatch.setenv("MULTIQC_MEGAQC_UPLOAD", "false")
    monkeypatch.setenv("MULTIQC_AI_SUMMARY", "false")
    monkeypatch.setenv("MULTIQC_AI_SUMMARY_FULL", "false")
