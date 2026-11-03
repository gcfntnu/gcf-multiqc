"""GCF summary-column compatibility over the upstream fastp parser."""

from multiqc.modules.fastp.fastp import MultiqcModule as NativeFastp


class MultiqcModule(NativeFastp):
    """Retain GCF's pre-filter read count and percent-passing-filter column IDs."""

    def find_log_files(self, sp_key, *args, **kwargs):
        if self.id == "gcf_fastp" and sp_key == "fastp":
            sp_key = "gcf_fastp"
        return super().find_log_files(sp_key, *args, **kwargs)

    def general_stats_addcols(self, data_by_sample, headers=None, *args, **kwargs):
        headers = dict(headers or {})
        # The legacy GCF report uses input reads as the main library-size
        # column. Native plots and raw exports still retain post-filter counts.
        headers.pop("filtering_result_passed_filter_reads", None)
        headers["before_filtering_total_reads"] = {
            "title": "Reads",
            "description": "Total reads before filtering",
            "scale": "GnBu",
            "shared_key": "read_count",
            "min": 0,
        }
        if "pct_surviving_reads" in headers:
            headers["pct_surviving"] = headers.pop("pct_surviving_reads")
        data_by_sample = {sample: dict(values) for sample, values in data_by_sample.items()}
        for values in data_by_sample.values():
            if "pct_surviving_reads" in values:
                values["pct_surviving"] = values["pct_surviving_reads"]
        return super().general_stats_addcols(data_by_sample, headers, *args, **kwargs)
