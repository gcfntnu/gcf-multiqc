"""GCF STAR columns with explicit, corrected splice-event semantics."""

from multiqc.modules.star.star import MultiqcModule as NativeSTAR
from multiqc.plots import bargraph


def add_gcf_metrics(data_by_sample):
    """Keep legacy keys while distinguishing mapped reads and splice events.

    Aligned reads follow upstream: unique plus multimapped, excluding reads
    mapping to too many loci. A read may contribute several splice events.
    Undefined ratios are absent values, not a misleading zero percentage.
    """
    enriched = {sample: dict(values) for sample, values in data_by_sample.items()}
    for values in enriched.values():
        values["total_mapped"] = values["mapped"]
        values["total_mapped_percent"] = values["mapped_percent"]
        splices = values.get("num_splices")
        unique = values.get("uniquely_mapped")
        annotated = values.get("num_annotated_splices")
        values["total_spliced_percent"] = 100 * splices / unique if splices is not None and unique else None
        values["novel_spliced_percent"] = (
            100 * (splices - annotated) / splices if splices and annotated is not None else None
        )
    return enriched


class MultiqcModule(NativeSTAR):
    """Preserve GCF column/section selectors using native parsing and plots.

    The historical ``total_spliced_percent`` key remains for configuration
    compatibility, but its label explicitly describes events per 100 unique
    reads. It is not the percentage of reads that are spliced and can exceed
    100. The former splice/read partition is replaced by an event-count plot.
    """

    def find_log_files(self, sp_key, *args, **kwargs):
        if self.id == "gcf_star" and (sp_key == "star" or sp_key.startswith("star/")):
            sp_key = "gcf_star" + sp_key.removeprefix("star")
        return super().find_log_files(sp_key, *args, **kwargs)

    def write_data_file(self, data, fn, *args, **kwargs):
        if fn == "multiqc_star":
            data = add_gcf_metrics(data)
        return super().write_data_file(data, fn, *args, **kwargs)

    def general_stats_addcols(self, data_by_sample, headers=None, *args, **kwargs):
        headers = dict(headers or {})
        headers["total_mapped"] = {
            **headers.pop("mapped"),
            "title": "Aligned reads",
            "hidden": False,
            "description": "Uniquely mapped plus multimapped reads; excludes too-many-loci reads",
        }
        headers["total_mapped_percent"] = {
            **headers.pop("mapped_percent"),
            "title": "% Aligned",
            "description": "Percentage of input reads uniquely mapped or multimapped; excludes too-many-loci reads",
        }
        headers["total_spliced_percent"] = {
            "title": "Splice events / 100 uniq reads",
            "description": "Splice events per 100 uniquely mapped reads; may exceed 100 because a read can span several junctions",
            "scale": "YlGn",
            "min": 0,
        }
        headers["novel_spliced_percent"] = {
            "title": "% Unannotated splice events",
            "description": "Percentage of splice events absent from the annotation; unavailable when there are no splice events",
            "suffix": "%",
            "scale": "YlGn",
            "min": 0,
            "max": 100,
        }
        return super().general_stats_addcols(add_gcf_metrics(data_by_sample), headers, *args, **kwargs)

    def stats_tables(self, data_by_sample):
        super().stats_tables(data_by_sample)
        self.add_section(
            name="Splice Events",
            anchor="star_spliced_reads",
            description="Counts of splice events by junction motif. A read can contribute multiple events.",
            helptext="These are splice events, not a partition of mapped reads. The aligned-read columns use upstream STAR semantics: uniquely mapped plus multimapped reads, excluding reads mapping to too many loci.",
            plot=bargraph.plot(
                data_by_sample,
                {
                    "num_GTAG_splices": {"name": "GT/AG"},
                    "num_GCAG_splices": {"name": "GC/AG"},
                    "num_ATAC_splices": {"name": "AT/AC"},
                    "num_noncanonical_splices": {"name": "Non-canonical"},
                },
                {
                    "id": "star_splice_plot",
                    "title": "STAR: Splice Events",
                    "ylab": "Splice events",
                    "cpswitch_counts_label": "Splice events",
                    "hide_zero_cats": False,
                },
            ),
        )
