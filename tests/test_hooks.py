"""Public selection/configuration contracts for the GCF compatibility hooks."""

from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from importlib_metadata import EntryPoint
from multiqc import config
import pytest
import yaml

from gcf_multiqc import hooks

FIXTURES = Path(__file__).parent / "fixtures"
ALIASES = {"gcf_fastp", "gcf_star", "gcf_cellranger_mkfastq"}


@pytest.fixture
def isolated_hook_config(monkeypatch):
    native = {
        name: EntryPoint(name=name, value=f"multiqc.modules.{name}.{name}:MultiqcModule", group="multiqc.modules.v1")
        for name in ("fastp", "star", "cellranger", "bclconvert")
    }
    native["custom_content"] = EntryPoint(
        name="custom_content",
        value="multiqc.core.special_case_modules.custom_content:custom_module_classes",
        group="multiqc.modules.v1",
    )
    monkeypatch.setattr(config, "avail_modules", dict(native))
    monkeypatch.setattr(
        config,
        "sp",
        {
            "kraken": {"fn": "*.kreport", "contents": "synthetic kraken"},
            "fastp": {"fn": "*.fastp.json"},
            "star": {"fn": "*Log.final.out"},
            "fastqc/data": {"fn": "fastqc_data.txt"},
        },
    )
    monkeypatch.setattr(config, "run_modules", [])
    monkeypatch.setattr(config, "exclude_modules", [])
    monkeypatch.setattr(config, "module_order", [])
    monkeypatch.setattr(config, "fn_clean_exts", [".fastq"])
    monkeypatch.setattr(config, "gcf_compatibility", True, raising=False)
    return native


def test_hook_preserves_user_search_rules_and_independent_alias_defaults(isolated_hook_config):
    custom = [{"fn": "my_parse_summary.csv", "contents": "site-specific signature"}]
    config.sp["parse"] = deepcopy(custom)
    hooks.before_config()
    assert config.sp["parse"] == custom
    assert config.sp["gcf_fastp"] == config.sp["fastp"]
    config.sp["gcf_fastp"]["fn"] = "alias-only.json"
    assert config.sp["fastp"]["fn"] == "*.fastp.json"
    assert config.sp["fastqc_rnaseq/data"] == config.sp["fastqc/data"]


def test_default_hooks_replace_legacy_names_without_automatic_alias_duplicates(isolated_hook_config):
    hooks.before_config()
    hooks.config_loaded()
    assert all(config.avail_modules[name].value.startswith("gcf_multiqc.") for name in isolated_hook_config)
    assert ALIASES <= set(config.exclude_modules)
    hooks.config_loaded()
    assert all(config.exclude_modules.count(alias) == 1 for alias in ALIASES)


def test_disabled_compatibility_preserves_native_entries_and_filename_rules(isolated_hook_config):
    config.gcf_compatibility = False
    config.run_modules = ["fastp", "star"]
    hooks.before_config()
    hooks.config_loaded()
    assert config.avail_modules == isolated_hook_config
    assert config.fn_clean_exts == [".fastq"]


@pytest.mark.parametrize("selection", ["cli", "module_order"])
def test_explicit_alias_remains_eligible_when_compatibility_is_disabled(isolated_hook_config, selection):
    config.gcf_compatibility = False
    if selection == "cli":
        config.run_modules = ["gcf_star"]
    else:
        config.module_order = [{"gcf_star": {"name": "Requested GCF STAR"}}]
    hooks.before_config()
    hooks.config_loaded()
    assert "gcf_star" not in config.exclude_modules
    assert config.avail_modules["star"] == isolated_hook_config["star"]


def _render(tmp_path, inputs, modules=(), overrides=None):
    settings = {
        "no_version_check": True,
        "megaqc_upload": False,
        "ai_summary": False,
        "ai_summary_full": False,
        "read_count_multiplier": 1,
        "read_count_prefix": "",
        **(overrides or {}),
    }
    config_file = tmp_path / "config.yaml"
    config_file.write_text(yaml.safe_dump(settings))
    command = [
        str(Path(sys.executable).with_name("multiqc")),
        str(inputs),
        "--strict",
        "--interactive",
        "--no-version-check",
        "--no-megaqc-upload",
        "--no-ai",
        "--config",
        str(config_file),
        "--outdir",
        str(tmp_path / "out"),
    ]
    if modules:
        command.append("--require-logs")
    for module in modules:
        command.extend(["-m", module])
    result = subprocess.run(
        command, cwd=tmp_path, env={**os.environ, "HOME": str(tmp_path)}, capture_output=True, text=True, timeout=60
    )
    (tmp_path / "render.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    return (
        json.loads((tmp_path / "out/multiqc_data/multiqc_data.json").read_text()),
        (tmp_path / "out/multiqc_report.html").read_text(),
    )


def _fastp_star_inputs(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for filename in ("Synthetic.fastp.json", "MultiSplice.Log.final.out"):
        shutil.copyfile(FIXTURES / "compatibility" / filename, inputs / filename)
    return inputs


def test_disabled_compatibility_cli_uses_native_fastp_and_star_metrics(tmp_path):
    data, html = _render(tmp_path, _fastp_star_inputs(tmp_path), ("fastp", "star"), {"gcf_compatibility": False})
    fastp_headers = data["report_general_stats_headers"]["fastp"]
    # Disabling compatibility also retains native naming, including the _R1 suffix.
    assert set(data["report_general_stats_data"]["fastp"]) == {"Synthetic_R1"}
    fastp = data["report_general_stats_data"]["fastp"]["Synthetic_R1"]
    assert fastp["filtering_result_passed_filter_reads"] == 800
    assert "filtering_result_passed_filter_reads" in fastp_headers
    assert "before_filtering_total_reads" not in fastp_headers
    assert "pct_surviving_reads" in fastp_headers
    assert "pct_surviving" not in fastp_headers
    star = data["report_saved_raw_data"]["multiqc_star"]["MultiSplice"]
    assert star["mapped"] == 70
    assert star["mapped_percent"] == 70
    assert "total_mapped" not in star
    assert "total_spliced_percent" not in star
    assert 'id="star_spliced_reads"' not in html


def test_default_cli_emits_legacy_mkfastq_once_without_explicit_alias(tmp_path):
    data, html = _render(tmp_path, FIXTURES / "singlecell/mkfastq")
    raw = data["report_saved_raw_data"]
    assert set(raw) == {"multiqc_cellranger"}
    assert set(raw["multiqc_cellranger"]) == {"001", "NA"}
    assert html.count('id="cellranger-mkfastq-qc"') == 1
    assert 'id="gcf_cellranger_mkfastq"' not in html


def test_default_cli_does_not_repeat_fastp_or_star_aliases(tmp_path):
    data, html = _render(tmp_path, _fastp_star_inputs(tmp_path))
    raw = data["report_saved_raw_data"]
    assert {name for name in raw if name.startswith("multiqc_fastp")} == {"multiqc_fastp"}
    assert {name for name in raw if name.startswith("multiqc_star")} == {"multiqc_star"}
    assert html.count('id="fastp-filtered-reads-chart"') == 1
    assert html.count('id="star_spliced_reads"') == 1
    assert data["report_general_stats_data"]["fastp"]["Synthetic"]["before_filtering_total_reads"] == 1000


def test_explicit_cli_aliases_keep_gcf_metrics_when_compatibility_is_disabled(tmp_path):
    data, html = _render(
        tmp_path, _fastp_star_inputs(tmp_path), ("gcf_fastp", "gcf_star"), {"gcf_compatibility": False}
    )
    assert data["report_general_stats_data"]["fastp"]["Synthetic_R1"]["before_filtering_total_reads"] == 1000
    assert "filtering_result_passed_filter_reads" not in data["report_general_stats_headers"]["fastp"]
    assert data["report_saved_raw_data"]["multiqc_star"]["MultiSplice"]["total_mapped"] == 70
    assert 'id="star_spliced_reads"' in html


def test_default_cli_honors_explicit_module_order_alias(tmp_path):
    data, html = _render(
        tmp_path,
        FIXTURES / "singlecell/mkfastq",
        overrides={
            "gcf_compatibility": False,
            "module_order": [{"gcf_cellranger_mkfastq": {}}],
            "exclude_modules": ["cellranger"],
        },
    )
    assert set(data["report_saved_raw_data"]) == {"multiqc_gcf_cellranger_mkfastq"}
    assert html.count('id="gcf_cellranger_mkfastq-mkfastq-qc"') == 1
