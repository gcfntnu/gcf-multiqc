# Synthetic native sequencing report inputs

Derived from the synthetic `RUN_INFO`, `INTEROP`, `write_convert` and
`write_bcl2fastq` fixtures in gcf-bfq `tests/test_sequencing_qc.py` at
`dca384260a53a874da08e66d78c5a73c6fb2a8dc`. The run identifier has been replaced
with an explicit synthetic value. No BFQ package, scientific tool, InterOp
executable, sequencer input, mail or container was used to create these files.

The BCL Convert directory deliberately colocates `RunInfo.xml` and
`Demultiplex_Stats.csv`, matching BFQ's isolated early-report staging layout.
There is one assigned cluster for s1, zero for s2 and 9,999 undetermined clusters.
Quality metrics are absent, so unavailable Q30 values must remain unavailable.

The bcl2fastq fixture has 75 assigned, zero and 25 undetermined clusters with
synthetic read-quality metrics. The InterOp fixture contains one summary and
four aggregate lane/read rows; a surface-specific row must not replace the
aggregate lane statistics. The index-summary placeholder from the BFQ tests is
not included because it is a mock string, not a valid native report.
