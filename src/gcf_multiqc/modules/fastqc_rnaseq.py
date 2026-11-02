"""GCF's FastQC selector, implemented by the upstream FastQC parser."""

from multiqc.modules.fastqc.fastqc import MultiqcModule as NativeFastQC
from multiqc.types import Anchor


class MultiqcModule(NativeFastQC):
    """Keep FastQC plots and the legacy selector without duplicate summary columns.

    The plugin registers ``fastqc_rnaseq/*`` search patterns before file search.
    The native parser remains responsible for FastQC data, assets and plots.
    """

    def __init__(self):
        self.mod_cust_config = {
            **self.mod_cust_config,
            "generalstats": False,
        }
        super().__init__()

    def find_log_files(self, sp_key, *args, **kwargs):
        # This is called immediately after native module initialization, before
        # parsing or plotting. Set our identity without inventing a custom
        # anchor, which would also rename the exported multiqc_fastqc file.
        self.anchor = Anchor(self.mod_cust_config.get("anchor", "fastqc_rnaseq"))
        self.intro = self._get_intro()
        if sp_key.startswith("fastqc/"):
            sp_key = "fastqc_rnaseq/" + sp_key.removeprefix("fastqc/")
        return super().find_log_files(sp_key, *args, **kwargs)

    def fastqc_general_stats(self):
        # Native FastQC may choose the Falco anchor from its input. The GCF
        # selector keeps a stable anchor for existing section-order settings.
        self.anchor = Anchor(self.mod_cust_config.get("anchor", "fastqc_rnaseq"))
        self.intro = self._get_intro()

    def add_section(self, *args, **kwargs):
        # User-supplied custom anchors are prepended by BaseMultiqcModule.
        # Otherwise use the stable GCF section IDs directly.
        prefix = "" if self.mod_cust_config.get("anchor") else "fastqc_rnaseq_"
        for key in ("id", "anchor"):
            value = kwargs.get(key)
            if isinstance(value, str):
                # Some native sections use self.anchor; others still use a
                # literal fastqc prefix. Normalize both exactly once.
                for native_prefix in (f"{self.anchor}_", "fastqc_"):
                    if value.startswith(native_prefix):
                        kwargs[key] = prefix + value.removeprefix(native_prefix)
                        break
        return super().add_section(*args, **kwargs)
