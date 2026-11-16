"""Small-RNA cutadapt reports (gcf-workflows smallrna, cutadapt filter route); not a general cutadapt module.

`gcf_smallrna_cutadapt` reads only the merged reports of the gcf-workflows smallrna cutadapt route
(gcf-workflows#350): one single-end JSON report per sample in cutadapt's schema, tagged "Cutadapt report, merged
from the passes in gcf_passes", with input reads from the 3' adapter pass, filtered categories summed over both
passes, and `--discard-trimmed` used for 5' adapter artefacts. Its search pattern matches that tag only, so per-pass
reports and text logs are never read. Reports select it explicitly (`-m gcf_smallrna_cutadapt`); selecting native
`cutadapt` as well is an error (hooks.config_loaded), and default runs leave it out. Paired-end reports and two
reports with one sample name raise an error.

Upstream MultiQC 1.35 reads only `too_short`, `too_long` and `too_many_n` from `read_counts.filtered` and shows the
rest as "Filtered reads (uncategorised)". This module names every JSON category, `discard_trimmed` as "5′ adapter
artefacts"; uncategorised remains for counts the report does not explain. General Statistics give the GCF fastp
module's columns (input reads, % passing filter, % with 3' adapter) instead of "Trimmed bases", which counts the bases
of discarded reads as trimmed.

Upstream coupling (MultiQC 1.35); compare with upstream's module before accepting a MultiQC bump:
- upstream's `__init__` finds files with `find_log_files("cutadapt")` (redirected here to this module's pattern) and
  writes General Statistics and the data file before it adds the plots;
- upstream's `parse_json` and private `_extract_sample_name_from_json` give the sample name and the `cutadapt_data`
  keys used here (`r_processed`, `r_written`, `r_too_short`, `r_too_long`, `r_too_many_N`, `r1_with_adapters`); each
  report is read a second time here for the other categories;
- upstream's `parse_file` recomputes `r_filtered_unexplained` for every sample after each file, so the recomputation
  here covers every sample too;
- `cutadapt_filtered_barplot` is re-implemented with upstream's category list, the extra categories and a small-RNA
  description.
"""

import json
import os

from multiqc.base_module import SampleGroupingConfig
from multiqc.modules.cutadapt.cutadapt import MultiqcModule as NativeCutadapt
from multiqc.plots import bargraph
from multiqc.types import ColumnKey, SampleName

MERGED_TAG = "Cutadapt report, merged from the passes in gcf_passes"
# JSON read_counts.filtered keys upstream does not read, with their names in the filtered-reads plot.
EXTRA_CATEGORIES = {
    "too_many_expected_errors": "Reads with too many expected errors",
    "casava_filtered": "Reads failing the CASAVA filter",
    "discard_trimmed": "5′ adapter artefacts",
    "discard_untrimmed": "Reads discarded: no adapter found",
}
UPSTREAM_CATEGORIES = ["r_too_short", "r_too_long", "r_too_many_N"]


class MultiqcModule(NativeCutadapt):
    """Merged small-RNA cutadapt reports: every JSON filtering category named; GCF General Statistics columns."""

    def __init__(self):
        self.gcf_sources = {}
        super().__init__()

    def find_log_files(self, sp_key, *args, **kwargs):
        if self.id == "gcf_smallrna_cutadapt" and sp_key == "cutadapt":
            sp_key = "gcf_smallrna_cutadapt"
        return super().find_log_files(sp_key, *args, **kwargs)

    def parse_log(self, f):
        raise ValueError(f"gcf_smallrna_cutadapt reads merged JSON reports only: {os.path.join(f['root'], f['fn'])}")

    def parse_json(self, f):
        source = os.path.join(f["root"], f["fn"])
        with open(source) as fh:
            data = json.load(fh)
        if data.get("tag") != MERGED_TAG or "gcf_passes" not in data:
            raise ValueError(f"Not a merged small-RNA cutadapt report: {source}")
        if data["input"]["paired"]:
            raise ValueError(f"Paired-end report {source}: the small-RNA merged reports are single-end")
        name = SampleName(self._extract_sample_name_from_json(data, f))
        if name in self.gcf_sources:
            raise ValueError(f"Duplicate small-RNA cutadapt sample {name!r}: {self.gcf_sources[name]} and {source}")
        self.gcf_sources[name] = source
        super().parse_json(f)
        d = self.cutadapt_data[name]
        for key in EXTRA_CATEGORIES:
            value = data["read_counts"]["filtered"].get(key)
            if value is not None:
                d["r_" + key] = value

    def parse_file(self, f):
        super().parse_file(f)
        # Upstream recomputes the remainder for every sample from its three categories; redo it with all of them.
        for d in self.cutadapt_data.values():
            if d.get("r_processed") is None:
                continue
            explained = sum(d.get(key) or 0 for key in [*UPSTREAM_CATEGORIES, "r_written"])
            explained += sum(d.get("r_" + key) or 0 for key in EXTRA_CATEGORIES)
            if d["r_processed"] - explained > 0:
                d["r_filtered_unexplained"] = d["r_processed"] - explained
            else:
                d.pop("r_filtered_unexplained", None)

    def cutadapt_general_stats_table(self):
        """Input reads, % passing filter and % with 3' adapter, as the GCF fastp columns; no "Trimmed bases"."""
        data = {}
        for s_name, d in self.cutadapt_data.items():
            values = {}
            if d.get("r_processed") is not None:
                values[ColumnKey("before_filtering_total_reads")] = d["r_processed"]
                if d["r_processed"] > 0:
                    if d.get("r_written") is not None:
                        values[ColumnKey("pct_surviving")] = 100.0 * d["r_written"] / d["r_processed"]
                    if d.get("r1_with_adapters") is not None:
                        values[ColumnKey("pct_adapter")] = 100.0 * d["r1_with_adapters"] / d["r_processed"]
            data[s_name] = values
        self.general_stats_addcols(
            data,
            {
                ColumnKey("before_filtering_total_reads"): {
                    "title": "Reads",
                    "description": "Total reads before filtering",
                    "scale": "GnBu",
                    "shared_key": "read_count",
                    "min": 0,
                },
                ColumnKey("pct_surviving"): {
                    "title": "% PF (Reads)",
                    "description": "Percent reads passing filter",
                    "suffix": "%",
                    "scale": "RdYlGn",
                    "min": 0,
                    "max": 100,
                },
                ColumnKey("pct_adapter"): {
                    "title": "% Adapter",
                    "description": "Percentage of reads in which the 3' adapter was found",
                    "suffix": "%",
                    "scale": "RdYlGn-rev",
                    "min": 0,
                    "max": 100,
                },
            },
            group_samples_config=SampleGroupingConfig(
                cols_to_weighted_average=[
                    (ColumnKey("pct_surviving"), ColumnKey("before_filtering_total_reads")),
                    (ColumnKey("pct_adapter"), ColumnKey("before_filtering_total_reads")),
                ],
            ),
        )

    def cutadapt_filtered_barplot(self):
        """Upstream's filtered-reads plot with every JSON category named."""
        cats = {
            "pairs_written": {"name": "Pairs passing filters"},
            "r_written": {"name": "Reads passing filters"},
            "pairs_too_short": {"name": "Pairs that were too short"},
            "r_too_short": {"name": "Reads that were too short"},
            "pairs_too_long": {"name": "Pairs that were too long"},
            "r_too_long": {"name": "Reads that were too long"},
            "pairs_too_many_N": {"name": "Pairs with too many N"},
            "r_too_many_N": {"name": "Reads with too many N"},
            **{"r_" + key: {"name": name} for key, name in EXTRA_CATEGORIES.items()},
            "pairs_filtered_unexplained": {"name": "Filtered pairs (uncategorised)"},
            "r_filtered_unexplained": {"name": "Filtered reads (uncategorised)"},
        }
        self.add_section(
            name="Filtered Reads",
            anchor="cutadapt_filtered_reads",
            description="This plot shows the number of reads removed by the small-RNA cutadapt filter, by reason.",
            plot=bargraph.plot(
                self.cutadapt_data,
                cats,
                {"id": "cutadapt_filtered_reads_plot", "title": "Cutadapt: Filtered Reads", "ylab": "Counts"},
            ),
        )
