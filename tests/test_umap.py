"""Render simultaneous legacy PCA/UMAP files through the installed plugin."""

import json
import os
from pathlib import Path
import subprocess
import sys

import yaml


def test_legacy_pca_and_umap_have_distinct_plots_and_literal_sample_metadata(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    identities = {
        "001": ("Synthetic Group A", "#ff0000"),
        "1": ("Synthetic Group B", "#0000ff"),
        "S1.fastq.gz": ("Synthetic Group A", "#ff0000"),
        "S1": ("Synthetic Group B", "#0000ff"),
    }
    originals = {}
    expected_points = {}
    for section, offset in (("pca", 0), ("umap", 10)):
        points = {
            sample: {"x": offset + index, "y": offset - index, "name": group, "color": color}
            for index, (sample, (group, color)) in enumerate(identities.items(), 1)
        }
        expected_points[section] = points
        content = {
            "id": section,
            "section_name": section.upper(),
            "parent_id": "exploratory",
            "parent_name": "Exploratory analysis",
            "parent_description": "Synthetic feature-table exploration.",
            "plot_type": "scatter",
            # The historical small-RNA producer emits these PCA settings for
            # UMAP too. Input correction belongs to the adapter, not the test.
            "pconfig": {
                "id": "pca_scatter",
                "title": "PCA",
                "xlab": "PC1",
                "ylab": "PC2",
                "tt_label": "PC1 {point.x:.2f}: PC2 {point.y:.2f}",
                "data_labels": [{"name": "PC1 vs PC2", "xlab": "PC1", "ylab": "PC2"}],
            },
            "data": [points],
        }
        fixture = inputs / f"{section}_mqc.yaml"
        fixture.write_text(yaml.safe_dump(content, sort_keys=False))
        originals[fixture] = fixture.read_bytes()

    settings = tmp_path / "config.yaml"
    settings.write_text("no_version_check: true\nmegaqc_upload: false\nai_summary: false\nai_summary_full: false\n")
    output = tmp_path / "output"
    command = [
        str(Path(sys.executable).with_name("multiqc")),
        str(inputs),
        "--strict",
        "--interactive",
        "--no-version-check",
        "--no-megaqc-upload",
        "--no-ai",
        "-m",
        "custom_content",
        "--config",
        str(settings),
        "--outdir",
        str(output),
    ]
    completed = subprocess.run(
        command, cwd=tmp_path, env={**os.environ, "HOME": str(tmp_path)}, capture_output=True, text=True, timeout=60
    )
    (tmp_path / "multiqc.log").write_text(completed.stdout + completed.stderr)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert all(fixture.read_bytes() == contents for fixture, contents in originals.items())
    assert (output / "multiqc_report.html").stat().st_size > 1000

    report = json.loads((output / "multiqc_data/multiqc_data.json").read_text())
    plots = report["report_plot_data"]
    assert {key for key in plots if key.endswith("_scatter")} == {"pca_scatter", "umap_scatter"}
    for section, prefix in (("pca", "PC"), ("umap", "UMAP")):
        plot = plots[f"{section}_scatter"]
        assert plot["pconfig"]["title"] == section.upper()
        assert plot["layout"]["showlegend"] is True
        assert len(plot["datasets"]) == 1
        dataset = plot["datasets"][0]
        assert plot["pconfig"]["xlab"] == f"{prefix}1"
        assert plot["pconfig"]["ylab"] == f"{prefix}2"
        assert dataset["layout"]["xaxis"]["title"]["text"] == f"{prefix}1"
        assert dataset["layout"]["yaxis"]["title"]["text"] == f"{prefix}2"
        points = {point["name"]: point for point in dataset["points"]}
        assert len(dataset["points"]) == len(identities)
        assert set(points) == set(identities)
        for sample, (group, color) in identities.items():
            assert points[sample]["group"] == group
            assert points[sample]["color"] == color
            assert points[sample]["x"] == expected_points[section][sample]["x"]
            assert points[sample]["y"] == expected_points[section][sample]["y"]
    umap = plots["umap_scatter"]["datasets"][0]
    assert umap["dconfig"]["title"] == "UMAP"
    assert "PC1" not in umap["trace_params"]["hovertemplate"]
    assert "PC2" not in umap["trace_params"]["hovertemplate"]
