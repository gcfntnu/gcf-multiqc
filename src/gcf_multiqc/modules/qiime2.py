"""GCF QIIME2 exports for MultiQC.

Adapted from gcfntnu/MultiQC's qiime2 module (GCF contributors, 2020).
Distributed under the GNU General Public License v3.0; see LICENSE.
"""

import colorsys
import csv
import hashlib
import io
import logging
import math
import re
from collections import defaultdict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, scatter

log = logging.getLogger(__name__)
TAXONOMIC_LEVELS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
DADA2_COLUMNS = ("input", "filtered", "denoised", "merged", "non-chimeric")


def _count(value):
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"Expected a finite nonnegative count, got {value!r}")
    return number


def parse_dada2(contents):
    """Return filtering losses by region and sample from a QIIME2 stats export."""
    reader = csv.DictReader(io.StringIO(contents), delimiter="\t")
    if not reader.fieldnames or not set(DADA2_COLUMNS).issubset(reader.fieldnames):
        raise ValueError("QIIME2 DADA2 metadata lacks required count columns")
    sample_column = reader.fieldnames[0]
    result = defaultdict(dict)
    for row in reader:
        sample_region = row[sample_column]
        if sample_region.startswith("#"):
            continue
        sample, separator, region = sample_region.rpartition("_")
        if not separator or not sample or not region:
            raise ValueError(f"Expected a sample_region identifier, got {sample_region!r}")
        cumulative = [_count(row[key]) for key in DADA2_COLUMNS]
        if any(after > before for before, after in zip(cumulative, cumulative[1:])):
            raise ValueError(f"DADA2 counts increase between filtering steps for {sample_region!r}")
        if sample in result[region]:
            raise ValueError(f"Duplicate DADA2 sample and region: {sample_region!r}")
        losses = {key: cumulative[i - 1] - cumulative[i] for i, key in enumerate(DADA2_COLUMNS) if i}
        result[region][sample] = {**losses, "passed": cumulative[-1]}
    return dict(result)


def parse_taxonomy(contents, filename):
    """Parse taxon columns only, retaining literal IDs and quoted CSV labels."""
    match = re.fullmatch(r"level-([1-7])\.csv", filename)
    if not match:
        raise ValueError(f"Unrecognized QIIME2 taxonomy filename: {filename!r}")
    level = TAXONOMIC_LEVELS[int(match[1]) - 1]
    reader = csv.reader(io.StringIO(contents))
    header = next(reader, [])
    taxa = []
    database = "custom"
    for index, label in enumerate(header[1:], 1):
        if label.startswith("D_0__"):
            database = "silva"
        elif label.startswith("k__Fungi"):
            database = "unite"
        elif label.startswith("k__Bacteria"):
            database = "greengenes"
        if re.match(r"(?:D_0__|[kd]__)", label):
            name = label.split(";")[-1].split("__", 1)[-1].strip() or "Unclassified"
            taxa.append((index, name))
    if not header or not taxa:
        raise ValueError(f"No taxonomy columns found in {filename!r}")
    data = {}
    for row in reader:
        if not row or row[0].startswith("#"):
            continue
        if len(row) != len(header) or not row[0]:
            raise ValueError(f"Malformed taxonomy row in {filename!r}")
        if row[0] in data:
            raise ValueError(f"Duplicate taxonomy sample: {row[0]!r}")
        counts = defaultdict(float)
        for index, name in taxa:
            counts[name] += _count(row[index])
        data[row[0]] = dict(counts)
    return level, data, database


def group_color(group):
    """Stable colors without a fixed palette or dependence on sample order."""
    digest = hashlib.sha256(group.encode("utf-8")).digest()
    hue = int.from_bytes(digest[:4], "big") / 2**32
    lightness = 0.42 + digest[4] / 255 * 0.16
    saturation = 0.55 + digest[5] / 255 * 0.3
    rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
    return "#" + "".join(f"{round(channel * 255):02x}" for channel in rgb)


def parse_ordination(contents):
    """Read bounded Site records and the workflow-appended Samples metadata."""
    lines = contents.splitlines()
    try:
        marker = next(i for i, line in enumerate(lines) if line.strip() == "Samples:")
    except StopIteration as exc:
        raise ValueError("QIIME2 ordination is missing the appended Samples: metadata marker") from exc
    site = next((i for i, line in enumerate(lines[:marker]) if line.startswith("Site\t")), None)
    if site is None:
        raise ValueError("QIIME2 ordination is missing the Site coordinates section")
    fields = lines[site].split("\t")
    if len(fields) != 3:
        raise ValueError("Malformed QIIME2 Site header")
    n_samples, n_axes = int(fields[1]), int(fields[2])
    if n_samples < 0 or n_axes < 2 or site + n_samples >= marker:
        raise ValueError("Truncated QIIME2 Site coordinates or fewer than two axes")
    points = {}
    for line in lines[site + 1 : site + 1 + n_samples]:
        values = line.split("\t")
        if len(values) != n_axes + 1 or not values[0]:
            raise ValueError("Malformed QIIME2 Site coordinate row")
        sample = values[0]
        if sample in points:
            raise ValueError(f"Duplicate ordination sample: {sample!r}")
        x, y = float(values[1]), float(values[2])
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError(f"Non-finite ordination coordinates for {sample!r}")
        points[sample] = {"x": x, "y": y, "name": sample}
    metadata = csv.DictReader(io.StringIO("\n".join(lines[marker + 1 :])), delimiter="\t")
    if not metadata.fieldnames or "Sample_ID" not in metadata.fieldnames:
        raise ValueError("QIIME2 appended metadata lacks the Sample_ID column")
    seen = set()
    for row in metadata:
        sample = row["Sample_ID"]
        if sample.startswith("#"):
            continue
        if sample in seen:
            raise ValueError(f"Duplicate ordination metadata sample: {sample!r}")
        seen.add(sample)
        group = row.get("Sample_Group", "")
        if sample in points and group:
            points[sample].update(group=group, color=group_color(group))
    return points


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="QIIME2", anchor="qiime2", href="https://qiime2.org/", info="GCF microbiome workflow summaries."
        )
        dada2 = defaultdict(dict)
        taxonomy = {}
        ordination = {}
        databases = set()
        for f in self.find_log_files("qiime2/dada2"):
            for region, samples in parse_dada2(f["f"]).items():
                for sample, values in samples.items():
                    if sample in dada2[region]:
                        raise ValueError(f"Duplicate DADA2 input for {sample!r}, region {region!r}")
                    dada2[region][sample] = values
                    self.add_data_source(f, sample, section="DADA2")
        for f in self.find_log_files("qiime2/taxonomy"):
            level, samples, database = parse_taxonomy(f["f"], f["fn"])
            if level == "Kingdom":
                continue  # Preserve the established report's Phylum–Species tabs.
            if level in taxonomy:
                raise ValueError(f"Duplicate QIIME2 taxonomy level: {level}")
            taxonomy[level] = self.ignore_samples(samples)
            databases.add(database)
            for sample in samples:
                self.add_data_source(f, sample, section="Taxonomy")
        for f in self.find_log_files("qiime2/rpca"):
            for sample, values in parse_ordination(f["f"]).items():
                if sample in ordination:
                    raise ValueError(f"Duplicate ordination input for {sample!r}")
                ordination[sample] = values
                self.add_data_source(f, sample, section="Ordination")
        dada2 = {region: self.ignore_samples(samples) for region, samples in dada2.items()}
        dada2 = {region: samples for region, samples in dada2.items() if samples}
        taxonomy = {level: samples for level, samples in taxonomy.items() if samples}
        ordination = self.ignore_samples(ordination)
        if not (dada2 or taxonomy or ordination):
            raise ModuleNoSamplesFound
        self.add_software_version(None)
        if dada2:
            total = defaultdict(lambda: defaultdict(float))
            for samples in dada2.values():
                for sample, values in samples.items():
                    for key, value in values.items():
                        total[sample][key] += value
            total = {sample: dict(values) for sample, values in total.items()}
            self.write_data_file({"Total": total, **dada2}, "multiqc_qiime2_dada2")
            self.general_stats_addcols(
                total,
                {
                    "passed": {
                        "title": "Dada2 PF PE reads",
                        "description": "Paired-end reads passing DADA2 filtering",
                        "min": 0,
                        "modify": lambda x: x * config.read_count_multiplier,
                        "scale": "GnBu",
                        "shared_key": "read_count",
                    }
                },
            )
            regions = sorted(dada2)
            datasets = [dada2[region] for region in regions]
            if len(regions) > 1:
                regions.insert(0, "Total")
                datasets.insert(0, total)
            categories = {
                key: {"name": label}
                for key, label in (
                    ("passed", "Passed Filter"),
                    ("filtered", "Low Quality"),
                    ("merged", "No Overlap"),
                    ("denoised", "Failed denoising"),
                    ("non-chimeric", "Chimeric"),
                )
            }
            self.add_section(
                name="Dada2 Filtered Reads",
                anchor="qiime2-dada2-stats",
                description="Filtering losses and retained paired-end reads for each amplified region.",
                plot=bargraph.plot(
                    datasets,
                    [categories] * len(datasets),
                    {
                        "id": "dada2-stats-plot",
                        "title": "Dada2: Sequence Stats",
                        "ylab": "# Reads",
                        "data_labels": regions,
                    },
                ),
            )
        if taxonomy:
            self.write_data_file(taxonomy, "multiqc_qiime2_taxonomy")
            levels = [level for level in TAXONOMIC_LEVELS if level in taxonomy]
            categories = []
            for level in levels:
                totals = defaultdict(float)
                for values in taxonomy[level].values():
                    for taxon, count in values.items():
                        totals[taxon] += count
                categories.append(
                    {taxon: {"name": taxon} for taxon in sorted(totals, key=lambda key: (-totals[key], key))}
                )
            self.add_section(
                name="Taxonomic Classification",
                anchor="qiime2-taxonomy",
                description=f"QIIME2 classifications using {', '.join(sorted(databases)).upper()} reference data.",
                plot=bargraph.plot(
                    [taxonomy[level] for level in levels],
                    categories,
                    {
                        "id": "taxonomy-classification",
                        "title": "QIIME2: Taxonomy",
                        "ylab": "# Reads",
                        "data_labels": levels,
                    },
                ),
            )
        if ordination:
            self.write_data_file(ordination, "multiqc_qiime2_rpca")
            self.add_section(
                name="Deicode",
                anchor="qiime2-deicode",
                description="DEICODE robust compositional PCA. Points retain sample IDs and are colored by Sample_Group when available.",
                plot=scatter.plot(
                    ordination,
                    {
                        "id": "dada2-rpca-plot",
                        "title": "QIIME2: Sample Ordination",
                        "xlab": "PC1",
                        "ylab": "PC2",
                        "showlegend": True,
                    },
                ),
            )
