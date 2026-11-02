"""Keep all local report tests offline, including their Python subprocesses."""

import os
from html.parser import HTMLParser
from pathlib import Path
import sys

import pytest

SUPPORT = Path(__file__).parent / "support"
sys.path.insert(0, str(SUPPORT))
import sitecustomize  # noqa: E402, F401


@pytest.fixture
def general_stats_cells():
    """Read visible text from one column of the actual rendered statistics table."""

    def read(html, column):
        class Cells(HTMLParser):
            active = False
            collecting = False
            sample = None

            def handle_starttag(self, tag, attributes):
                attributes = dict(attributes)
                if tag == "table":
                    self.active = attributes.get("id") == "general_stats_table_table"
                if self.active:
                    if tag == "tr":
                        self.sample = attributes.get("data-sample-group")
                    elif tag == "td" and column in attributes.get("class", "").split():
                        self.collecting = True
                        values[self.sample] = ""

            def handle_data(self, text):
                if self.active and self.collecting:
                    values[self.sample] += text

            def handle_endtag(self, tag):
                if tag == "table":
                    self.active = False
                elif tag == "td":
                    self.collecting = False

        values = {}
        Cells().feed(html)
        assert values, f"No rendered General Statistics cells for {column}"
        return values

    return read


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
