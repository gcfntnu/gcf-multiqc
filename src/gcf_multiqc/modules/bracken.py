"""Restore the GCF Bracken selector using MultiQC's current Kraken parser."""

from multiqc.modules.kraken.kraken import MultiqcModule as NativeKraken


class MultiqcModule(NativeKraken):
    """Report staged Bracken kreports separately from their Kraken inputs.

    The plugin's search patterns identify Bracken by its report filename and
    exclude those files from Kraken when both modules run. A missing
    unclassified row is not sufficient to identify a Bracken report: valid
    Kraken reports can omit it too. The seven-column Bracken abundance table
    is a different format and is not parsed by this module.
    """

    def __init__(self):
        super().__init__(
            name="Bracken",
            anchor="bracken",
            href="https://ccb.jhu.edu/software/bracken/",
            info="Estimates taxonomic abundance from Kraken classifications.",
            doi="10.7717/peerj-cs.104",
            sp_key="bracken",
        )
