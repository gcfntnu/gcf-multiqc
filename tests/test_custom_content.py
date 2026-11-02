import json
from pathlib import Path
import subprocess
import sys

import yaml


def render(tmp_path, content, config=None):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fixture = inputs / "fixture_mqc.yaml"
    original = yaml.safe_dump(content)
    fixture.write_text(original)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config or {}))
    result = subprocess.run([str(Path(sys.executable).with_name("multiqc")), "--strict", "--no-version-check",
                             "--no-megaqc-upload", "--no-ai", "--interactive", "-m", "custom_content",
                             "-c", str(config_path), "-o", str(tmp_path / "out"), str(inputs)],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert fixture.read_text() == original
    assert (tmp_path / "out/multiqc_report.html").stat().st_size > 1000
    return json.loads((tmp_path / "out/multiqc_data/multiqc_data.json").read_text())


def test_custom_literal_ids_do_not_collapse(tmp_path):
    data = {"S1.fastq.gz": {"value": 1}, "S1": {"value": 2}, "001": {"value": 3}, "1": {"value": 4}}
    report = render(tmp_path, {"id": "literal", "section_name": "Literal IDs", "plot_type": "table", "data": data})
    assert report["report_saved_raw_data"]["multiqc_literal"] == data


def test_pca_sample_identity_and_dataset_labels(tmp_path):
    datasets = [{"001": {"x": 1, "y": 2, "name": "Group A", "color": "#ff0000"}},
                {"001": {"x": 3, "y": 4, "name": "Group A", "color": "#ff0000"}}]
    report = render(tmp_path, {"id": "pca", "section_name": "PCA", "plot_type": "scatter", "data": datasets,
                              "pconfig": {"id": "pca_scatter", "title": "PCA", "tt_label": "PC1 {point.x}: PC2 {point.y}",
                                          "data_labels": [{"name": "PC1 vs PC2", "xlab": "PC1", "ylab": "PC2"},
                                                          {"name": "PC3 vs PC4", "xlab": "PC3", "ylab": "PC4", "xmax": 5}]}})
    plot = report["report_plot_data"]["pca_scatter"]
    assert len(plot["datasets"]) == 2
    for dataset in plot["datasets"]:
        assert dataset["points"][0]["name"] == "001"
        assert dataset["points"][0]["group"] == "Group A"
        assert dataset["points"][0]["color"] == "#ff0000"
    assert plot["datasets"][1]["dconfig"]["xlab"] == "PC3"
    assert plot["datasets"][1]["dconfig"]["xmax"] == 5
    assert "PC1" not in plot["datasets"][1]["trace_params"]["hovertemplate"]


def test_generalstats_preserves_text_metadata(tmp_path):
    data = {"001": {"External_ID": "0007", "Concentration": "<1"}, "1": {"External_ID": "NA", "Concentration": 12.35}}
    config = {"custom_data": {"general_statistics": {"plot_type": "generalstats", "pconfig": [{
        "External_ID": {"namespace": "gcf", "format": "{}"}, "Concentration": {"namespace": "gcf", "format": "{}"}}], "data": data}}}
    report = render(tmp_path, {"id": "summary", "section_name": "Summary", "plot_type": "html", "data": "<p>Synthetic</p>"}, config)
    assert report["report_general_stats_data"]["custom_content"] == data
