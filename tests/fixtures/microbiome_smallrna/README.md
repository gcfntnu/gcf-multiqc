# Synthetic microbiome and small RNA reports

All identifiers, counts and labels in this directory are invented. These files
model the layouts collected by GCF workflows; no scientific tool was executed.
The identifiers `001`, `1` and `NA` deliberately exercise literal sample names.

QIIME2 metadata records cumulative DADA2 counts and sample-region suffixes.
Taxonomy CSVs include metadata columns and quoted taxonomy labels. The ordination
contains the `Samples:` table appended by the GCF workflow.

Unitas annotation children are indented subsets, not extra reads. Its samples
have unequal totals, and the miRNA table includes a gene exactly at the detection
threshold. Histogram files retain the producer's two-line preamble.
