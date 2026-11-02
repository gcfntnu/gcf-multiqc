"""Register compatibility through MultiQC's supported hooks, without changing its files."""

from copy import deepcopy
from importlib.resources import files
import logging

from importlib_metadata import EntryPoint
from multiqc import config
import yaml

log = logging.getLogger(__name__)


def before_config():
    """Defaults are installed before explicit user config, which remains authoritative."""
    patterns = yaml.safe_load(files("gcf_multiqc").joinpath("search_patterns.yaml").read_text())
    for key, value in patterns.items():
        config.sp.setdefault(key, value)
    # Bracken uses the Kraken report format, but needs an independent identity.
    for suffix in ["*.bracken", "*.bracken_kreport", "*.bracken.kreport"]:
        pattern = deepcopy(config.sp["kraken"])
        pattern["fn"] = suffix
        config.sp.setdefault("bracken", []).append(pattern)
    for source, alias in [("fastqc", "fastqc_rnaseq"), ("fastp", "gcf_fastp"), ("star", "gcf_star")]:
        for key, value in list(config.sp.items()):
            if key == source or key.startswith(source + "/"):
                config.sp.setdefault(alias + key[len(source) :], deepcopy(value))


def config_loaded():
    """Preserve legacy selections deterministically; no duplicate entry point names."""
    # Prevent Kraken consuming Bracken files first (or reporting them twice).
    # Limit this partition to reports that actually enable the Bracken adapter.
    if (not config.run_modules or "bracken" in config.run_modules) and "bracken" not in config.exclude_modules:
        patterns = config.sp["kraken"]
        for pattern in patterns if isinstance(patterns, list) else [patterns]:
            excluded = pattern.setdefault("exclude_fn", [])
            if isinstance(excluded, str):
                excluded = pattern["exclude_fn"] = [excluded]
            for bracken_pattern in config.sp.get("bracken", []):
                filename = bracken_pattern.get("fn")
                if filename and filename not in excluded:
                    excluded.append(filename)
    if not getattr(config, "gcf_compatibility", True):
        return
    replacements = {
        "cellranger": "gcf_multiqc.modules.cellranger_mkfastq:MultiqcModule",
        "fastp": "gcf_multiqc.modules.fastp:MultiqcModule",
        "star": "gcf_multiqc.modules.star:MultiqcModule",
        "custom_content": "gcf_multiqc.custom_content:custom_module_classes",
    }
    for name, target in replacements.items():
        config.avail_modules[name] = EntryPoint(name=name, value=target, group="multiqc.modules.v1")
    # These are the historical fork's defaults. Existing workflow configuration
    # adds its own cleaning rules; the custom-content adapter preserves literal IDs.
    for extension in [".namesorted", "_R1", "_R2"]:
        if extension not in config.fn_clean_exts:
            config.fn_clean_exts.append(extension)
    log.info("GCF compatibility enabled (gcf-multiqc); upstream engine %s", config.version)
