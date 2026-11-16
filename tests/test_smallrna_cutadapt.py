"""gcf_smallrna_cutadapt: merged small-RNA cutadapt reports (gcf-workflows smallrna cutadapt route) with every JSON
filtering category named and GCF General Statistics columns (gcfntnu/gcf-multiqc#9). Native `cutadapt` stays
upstream's.

Fixtures are cutadapt 5.2 reports of synthetic reads. `Synthetic.cutadapt.json` merges the two passes of the
small-RNA route (3' adapter; then 5' adapter fragments discarded, minimum length, N limit, expected errors).
`Paired.cutadapt.json` (a paired-end pass) and `Synthetic.cutadapt.log` (the filter pass's console report) are
reports the module must not read.
"""

import json
from pathlib import Path
import shutil

import pytest

from test_compatibility import render, rendered_columns

FIXTURES = Path(__file__).parent / "fixtures" / "cutadapt"
MODULE = ["gcf_smallrna_cutadapt"]

CATEGORIES = [
    "r_written",
    "r_too_short",
    "r_too_long",
    "r_too_many_N",
    "r_too_many_expected_errors",
    "r_casava_filtered",
    "r_discard_trimmed",
    "r_discard_untrimmed",
    "r_filtered_unexplained",
]


def category_names(data):
    return [cat["name"] for cat in data["report_plot_data"]["cutadapt_filtered_reads_plot"]["datasets"][0]["cats"]]


def merged_report():
    return json.loads((FIXTURES / "Synthetic.cutadapt.json").read_text())


def write(folder, name, payload):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / name).write_text(json.dumps(payload, indent=2))


def inputs(tmp_path, *names):
    folder = tmp_path / "inputs"
    folder.mkdir(exist_ok=True)
    for name in names:
        shutil.copy(FIXTURES / name, folder)
    return folder


def test_categories_are_named_and_add_up_to_the_input(tmp_path):
    data, _, _ = render(tmp_path, inputs(tmp_path, "Synthetic.cutadapt.json"), MODULE)
    raw = data["report_saved_raw_data"]["multiqc_cutadapt"]["Synthetic"]
    assert raw["r_processed"] == 1000 and raw["r_written"] == 700
    assert raw["r_too_short"] == 150 and raw["r_too_many_N"] == 30
    assert raw["r_too_many_expected_errors"] == 20 and raw["r_discard_trimmed"] == 100
    assert "r_filtered_unexplained" not in raw
    assert sum(raw.get(key) or 0 for key in CATEGORIES) == raw["r_processed"]
    assert category_names(data) == [
        "Reads passing filters",
        "Reads that were too short",
        "Reads with too many N",
        "Reads with too many expected errors",
        "5′ adapter artefacts",
    ]


def test_every_sample_keeps_all_its_categories(tmp_path):
    """Upstream recomputes the uncategorised remainder for every sample after each file."""
    second = merged_report()
    second["input"]["path1"] = "Second.fastq"
    second["read_counts"]["filtered"].update(discard_trimmed=50, casava_filtered=10)
    second["read_counts"]["output"] = 740
    write(tmp_path / "inputs" / "a", "Synthetic.cutadapt.json", merged_report())
    write(tmp_path / "inputs" / "b", "Second.cutadapt.json", second)
    data, _, _ = render(tmp_path, tmp_path / "inputs", MODULE)
    raw = data["report_saved_raw_data"]["multiqc_cutadapt"]
    assert sorted(raw) == ["Second", "Synthetic"]
    for values in raw.values():
        assert "r_filtered_unexplained" not in values
        assert sum(values.get(key) or 0 for key in CATEGORIES) == values["r_processed"] == 1000
    assert raw["Second"]["r_casava_filtered"] == 10 and raw["Second"]["r_discard_trimmed"] == 50
    assert "Reads failing the CASAVA filter" in category_names(data)


def test_general_statistics_are_the_gcf_filter_columns(tmp_path):
    data, _, _ = render(
        tmp_path,
        inputs(tmp_path, "Synthetic.cutadapt.json"),
        MODULE,
        {"table_columns_placement": {"cutadapt": {"before_filtering_total_reads": 110, "pct_surviving": 112}}},
    )
    values = data["report_general_stats_data"]["cutadapt"]["Synthetic"]
    assert values == {"before_filtering_total_reads": 1000, "pct_surviving": 70, "pct_adapter": 90}
    headers = data["report_general_stats_headers"]["cutadapt"]
    assert "percent_trimmed" not in headers
    assert headers["before_filtering_total_reads"]["shared_key"] == "read_count"
    columns = rendered_columns(data, "cutadapt")
    assert columns["before_filtering_total_reads"]["placement"] == 110
    assert columns["pct_surviving"]["placement"] == 112


def test_only_merged_reports_are_read(tmp_path):
    """Per-pass reports and logs of the same samples, as in the workflow's provenance folder, are not read."""
    folder = inputs(tmp_path, "Synthetic.cutadapt.json", "Synthetic.cutadapt.log", "Paired.cutadapt.json")
    single_pass = merged_report()
    single_pass["tag"] = "Cutadapt report"
    del single_pass["gcf_passes"]
    write(folder / "provenance", "Synthetic.filter.json", single_pass)
    data, _, _ = render(tmp_path, folder, MODULE)
    raw = data["report_saved_raw_data"]["multiqc_cutadapt"]
    assert list(raw) == ["Synthetic"] and raw["Synthetic"]["r1_with_adapters"] == 900
    assert data["report_general_stats_data"]["cutadapt"]["Synthetic"]["pct_adapter"] == 90


def test_paired_end_report_fails(tmp_path):
    paired = merged_report()
    paired["input"].update(paired=True, path2="Synthetic_R2.fastq")
    write(tmp_path / "inputs", "Synthetic.cutadapt.json", paired)
    with pytest.raises(AssertionError, match="Paired-end report"):
        render(tmp_path, tmp_path / "inputs", MODULE)


def test_duplicate_sample_names_fail(tmp_path):
    write(tmp_path / "inputs" / "a", "Synthetic.cutadapt.json", merged_report())
    write(tmp_path / "inputs" / "b", "Synthetic.cutadapt.json", merged_report())
    with pytest.raises(AssertionError, match="Duplicate small-RNA cutadapt sample 'Synthetic'"):
        render(tmp_path, tmp_path / "inputs", MODULE)


def test_selecting_native_cutadapt_too_fails(tmp_path):
    with pytest.raises(AssertionError, match="Select cutadapt or gcf_smallrna_cutadapt, not both"):
        render(tmp_path, inputs(tmp_path, "Synthetic.cutadapt.json"), ["cutadapt", *MODULE])


@pytest.mark.parametrize("compatibility", [True, False])
def test_native_cutadapt_is_unchanged(tmp_path, compatibility):
    """Only explicit gcf_smallrna_cutadapt selections use this module; `cutadapt` stays upstream's."""
    data, _, _ = render(
        tmp_path, inputs(tmp_path, "Synthetic.cutadapt.json"), ["cutadapt"], {"gcf_compatibility": compatibility}
    )
    raw = data["report_saved_raw_data"]["multiqc_cutadapt"]["Synthetic"]
    assert raw["r_filtered_unexplained"] == 120 and "r_discard_trimmed" not in raw
    assert "percent_trimmed" in data["report_general_stats_headers"]["cutadapt"]
    assert "Filtered reads (uncategorised)" in category_names(data)


def test_zero_input_has_no_invented_percentages(tmp_path):
    payload = merged_report()
    counts = payload["read_counts"]
    counts.update(input=0, output=0, read1_with_adapter=0)
    counts["filtered"] = {key: (0 if value is not None else None) for key, value in counts["filtered"].items()}
    payload["basepair_counts"].update(input=0, input_read1=0, output=0, output_read1=0)
    write(tmp_path / "inputs", "Synthetic.cutadapt.json", payload)
    data, _, _ = render(tmp_path, tmp_path / "inputs", MODULE)
    assert data["report_general_stats_data"]["cutadapt"]["Synthetic"] == {"before_filtering_total_reads": 0}
