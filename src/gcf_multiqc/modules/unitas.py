"""GCF Unitas summaries for MultiQC.

Adapted from gcfntnu/MultiQC's unitas module (GCF contributors, 2020–2021).
Distributed under the GNU General Public License v3.0; see LICENSE.
"""

import logging
import math
from collections import defaultdict
from pathlib import Path
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)
BIOTYPES = ("miRNA", "tRNA", "rRNA", "protein_coding", "lincRNA", "snoRNA", "no annotation")


def _count(value):
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"Expected a finite nonnegative count, got {value!r}")
    return number


def parse_annotation(contents):
    """Count top-level annotations once; indented children are subsets."""
    annotations = defaultdict(float)
    total = 0.0
    for line in contents.splitlines():
        if not line.strip() or line.startswith("#") or line[0].isspace():
            continue
        fields = line.rsplit("\t", 1)
        if len(fields) != 2:
            raise ValueError("Malformed Unitas annotation row; expected category and count")
        label, raw_count = fields
        count = _count(raw_count)
        category = label.strip() if label.strip() in BIOTYPES else "other"
        annotations[category] += count
        total += count
    if not annotations:
        return {}, {}
    metrics = {"miRNA Percentage": annotations.get("miRNA", 0) / total if total else None}
    return dict(annotations), metrics


def sequence_biotype(filename):
    if re.fullmatch(r"unitas\.miR\..+\.info", filename):
        return "miRNA"
    for biotype in BIOTYPES:
        if filename == f"unitas.{biotype.replace(' ', '-')}.info":
            return biotype
    return None


def parse_sequence_lengths(contents):
    """Read the two-header-line Unitas length histogram without zero filling."""
    result = defaultdict(float)
    for line in contents.splitlines()[2:]:
        if not line.strip():
            break
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError("Malformed Unitas length histogram row")
        length = int(fields[0])
        if length < 1:
            raise ValueError("Unitas sequence length must be positive")
        result[length] += _count(fields[1] or "0")
    return dict(result)


def parse_mirna_counts(contents, min_detect=3):
    """Count miRNAs with at least three reads, allowing the producer's headers."""
    counts = {}
    header_lines = 0
    for line in contents.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        try:
            if len(fields) != 2:
                raise ValueError("Expected miRNA identifier and count")
            float(fields[1])
        except ValueError:
            if not counts and header_lines < 2:
                header_lines += 1
                continue
            raise ValueError("Malformed Unitas miRNA count row") from None
        count = _count(fields[1])
        if fields[0] in counts:
            raise ValueError(f"Duplicate Unitas miRNA: {fields[0]!r}")
        counts[fields[0]] = count
    return sum(value >= min_detect for value in counts.values()) if counts else None


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="unitas",
            anchor="unitas",
            href="https://www.smallrnagroup.uni-mainz.de/software.html",
            info="Annotation and length summaries from the GCF small-RNA workflow.",
        )
        annotations = {}
        metrics = {}
        lengths = defaultdict(dict)
        for f in self.find_log_files("unitas/annotation"):
            sample = Path(f["root"]).name
            sample_annotations, sample_metrics = parse_annotation(f["f"])
            if sample_annotations:
                if sample in annotations:
                    raise ValueError(f"Duplicate Unitas annotation sample: {sample!r}")
                annotations[sample] = sample_annotations
                metrics[sample] = sample_metrics
                self.add_data_source(f, sample, section="Annotations")
        for f in self.find_log_files("unitas/seqlen"):
            biotype = sequence_biotype(f["fn"])
            if biotype is None:
                continue
            sample = Path(f["root"]).parent.name
            histogram = parse_sequence_lengths(f["f"])
            if histogram:
                existing = lengths[biotype].setdefault(sample, {})
                for length, count in histogram.items():
                    existing[length] = existing.get(length, 0) + count
                self.add_data_source(f, sample, section="Sequence lengths")
        mirna_samples = set()
        for f in self.find_log_files("unitas/mirna"):
            sample = Path(f["root"]).name
            detected = parse_mirna_counts(f["f"])
            if detected is not None:
                if sample in mirna_samples:
                    raise ValueError(f"Duplicate Unitas miRNA count sample: {sample!r}")
                mirna_samples.add(sample)
                metrics.setdefault(sample, {})["miRNA Detected"] = detected
                self.add_data_source(f, sample, section="miRNA counts")
        annotations = self.ignore_samples(annotations)
        metrics = self.ignore_samples(metrics)
        lengths = {biotype: self.ignore_samples(samples) for biotype, samples in lengths.items()}
        lengths = {biotype: samples for biotype, samples in lengths.items() if samples}
        if not (annotations or metrics or lengths):
            raise ModuleNoSamplesFound
        self.add_software_version(None)
        if metrics:
            self.write_data_file(metrics, "multiqc_unitas")
            self.general_stats_addcols(
                metrics,
                {
                    "miRNA Percentage": {
                        "title": "% miRNA",
                        "description": "Fraction of top-level annotated reads classified as miRNA; unavailable for zero reads",
                        "max": 100,
                        "min": 0,
                        "suffix": "%",
                        "scale": "YlGn",
                        "modify": lambda x: x * 100,
                    },
                    "miRNA Detected": {
                        "title": "# Genes",
                        "description": "miRNA genes detected with at least 3 reads",
                        "min": 0,
                        "scale": "Blues",
                        "format": "{:,.0f}",
                    },
                },
            )
        if annotations:
            self.write_data_file(annotations, "multiqc_unitas_annotations")
            categories = {
                key: {"name": key} for key in (*BIOTYPES, "other") if any(key in row for row in annotations.values())
            }
            self.add_section(
                name="Annotations",
                anchor="annotations",
                description="Top-level annotation categories; subcategories are not counted twice.",
                plot=bargraph.plot(
                    annotations,
                    categories,
                    {
                        "id": "unitas_annotation",
                        "title": "Unitas: Gene annotations",
                        "ylab": "# Reads",
                        "cpswitch_c_active": True,
                    },
                ),
            )
        if lengths:
            self.write_data_file(lengths, "multiqc_unitas_lengths")
            biotypes = [key for key in BIOTYPES if key in lengths]
            self.add_section(
                name="Sequence Length Distribution",
                anchor="unitas_seqlen",
                plot=linegraph.plot(
                    [lengths[key] for key in biotypes],
                    {
                        "id": "unitas_seqlen_plot",
                        "title": "Unitas: Sequence length distribution",
                        "ylab": "# Reads",
                        "xlab": "Sequence length (nt)",
                        "xmin": 10,
                        "xmax": 50,
                        "data_labels": [{"name": key, "ylab": f"Number of {key} Reads"} for key in biotypes],
                    },
                ),
            )
