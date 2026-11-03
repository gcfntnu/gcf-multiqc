"""Preserve literal custom-content sample IDs while using the upstream renderer."""

from copy import deepcopy
import json
from pathlib import Path

from multiqc import config, report
from multiqc.base_module import BaseMultiqcModule
from multiqc.core.special_case_modules.custom_content import custom_module_classes as upstream_custom_content
from multiqc.plots.table_object import Cell
import yaml


def adapt_scatter(content):
    """Translate the existing GCF PCA/UMAP convention without recalculating data."""
    if (
        not isinstance(content, dict)
        or content.get("plot_type") != "scatter"
        or content.get("id") not in {"pca", "umap"}
    ):
        return None
    content = deepcopy(content)
    data = content.get("data", {})
    datasets = data if isinstance(data, list) else [data]
    for dataset in datasets:
        if not isinstance(dataset, dict):
            raise ValueError("GCF PCA/UMAP data must map sample IDs to points")
        for sample, points in dataset.items():
            for point in points if isinstance(points, list) else [points]:
                if isinstance(point, dict):
                    # Workflow producers use name for the GROUP, not the sample.
                    if point.get("name") is not None and str(point["name"]) != str(sample):
                        point.setdefault("group", str(point["name"]))
                    point["name"] = str(sample)
    pconfig = content.setdefault("pconfig", {})
    pconfig.setdefault("showlegend", True)
    # The old small-RNA producer reuses pca_scatter for UMAP as well.
    if content["id"] == "umap" and pconfig.get("id") == "pca_scatter":
        pconfig.update(id="umap_scatter", title="UMAP", xlab="UMAP1", ylab="UMAP2")
        for index, label in enumerate(pconfig.get("data_labels", [])):
            if isinstance(label, dict):
                label.update(
                    name=f"UMAP{2 * index + 1} vs UMAP{2 * index + 2}",
                    xlab=f"UMAP{2 * index + 1}",
                    ylab=f"UMAP{2 * index + 2}",
                )
    # Do not copy the fixed PC1/PC2 tooltip onto every component pair.
    pconfig.pop("tt_label", None)
    return content


def _stage_scatter_inputs():
    """Adapt recognized files in memory; never rewrite a workflow's input files."""
    scanner = BaseMultiqcModule(name="Custom content", anchor="custom_content")
    selected = []
    for item in scanner.find_log_files("custom_content"):
        suffix = Path(item["fn"]).suffix.lower()
        if suffix not in {".yaml", ".yml", ".json"}:
            continue
        try:
            content = json.loads(item["f"]) if suffix == ".json" else yaml.safe_load(item["f"])
        except (ValueError, yaml.YAMLError):
            continue  # Leave diagnostics for unrecognized content to upstream.
        adapted = adapt_scatter(content)
        if adapted is None:
            continue
        section = adapted["id"]
        if section in config.custom_data and "data" in config.custom_data[section]:
            raise ValueError(f"Multiple GCF scatter inputs define section {section!r}; use distinct section IDs")
        overrides = config.custom_data.get(section, {})
        config.custom_data[section] = {**adapted, **overrides}
        selected.append((item["root"], item["fn"]))
        for dataset in adapted["data"] if isinstance(adapted["data"], list) else [adapted["data"]]:
            for sample in dataset:
                scanner.add_data_source(item, s_name=str(sample), section=section)
    if selected:
        report.files["custom_content"] = [
            item for item in report.files.get("custom_content", []) if (item["root"], item["fn"]) not in selected
        ]


def _preserve_general_stats_text(sections):
    """Honor GCF's plain-text format without disabling numeric metric parsing."""
    for section in sections:
        columns = {
            key
            for key, header in report.general_stats_headers[section].items()
            if header.get("format") == "{}" and not header.get("modify")
        }
        for rows in report.general_stats_data[section].values():
            for row in rows:
                for key in columns & row.data.keys():
                    value = row.data[key]
                    if isinstance(value, str) and value.strip():
                        # General Statistics reparses strings as numbers even with
                        # format="{}". Upstream Cell bypasses that coercion and its
                        # report exporter preserves Cell.raw as the original value.
                        # Character references also protect punctuation from the
                        # renderer's unconditional numeric locale substitutions.
                        literal_html = "".join(f"&#{ord(character)};" for character in value)
                        row.data[key] = Cell(raw=value, mod=value, fmt=literal_html)


def custom_module_classes():
    # GCF's historical custom-content reader did not clean sample-key dictionaries.
    # Scope this compatibility behavior to custom content; standard tool parsers
    # still need filename cleaning to join their measurements to the same sample.
    clean_names = config.fn_clean_sample_names
    custom_data = config.custom_data
    custom_files = report.files.get("custom_content")
    config.fn_clean_sample_names = False
    config.custom_data = deepcopy(custom_data)
    previous_sections = set(report.general_stats_data)
    try:
        _stage_scatter_inputs()
        modules = upstream_custom_content()
        _preserve_general_stats_text(report.general_stats_data.keys() - previous_sections)
        return modules
    finally:
        config.fn_clean_sample_names = clean_names
        config.custom_data = custom_data
        if custom_files is not None:
            report.files["custom_content"] = custom_files
