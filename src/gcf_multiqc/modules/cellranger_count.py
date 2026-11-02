"""Cell Ranger count CSV adapter, adapted from gcfntnu/MultiQC (GPL-3.0).

Metric keys and units are retained. CSV count/percentage strings are normalized
to numbers in exports as well as plots; percentage fields are percentage points,
not fractions. Literal sample identifiers come only from filenames/directories.
"""

from __future__ import annotations

import csv
import logging
import math
from pathlib import Path
from typing import TextIO

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)

COUNTS = {
    "Number of Reads",
    "Estimated Number of Cells",
    "Mean Reads per Cell",
    "Median UMI Counts per Cell",
    "Median Genes per Cell",
    "Total Genes Detected",
}
PERCENTAGES = {
    "Valid Barcodes",
    "Sequencing Saturation",
    "Q30 Bases in Barcode",
    "Q30 Bases in RNA Read",
    "Reads Mapped to Genome",
    "Reads Mapped Confidently to Genome",
    "Reads Mapped Confidently to Transcriptome",
    "Reads Mapped Confidently to Exonic Regions",
    "Reads Mapped Confidently to Intronic Regions",
    "Reads Mapped Confidently to Intergenic Regions",
}
REGIONS = {
    "Reads Mapped Confidently to Exonic Regions": {"name": "Exonic"},
    "Reads Mapped Confidently to Intronic Regions": {"name": "Intronic"},
    "Reads Mapped Confidently to Intergenic Regions": {"name": "Intergenic"},
}


def parse_summary(handle: TextIO) -> dict:
    rows = [row for row in csv.reader(handle) if row and any(value.strip() for value in row)]
    if len(rows) != 2:
        raise ValueError("Cell Ranger count CSV must contain exactly one header and one data row")
    header, values = rows
    header[0] = header[0].removeprefix("\ufeff")
    if len(header) != len(values) or not header:
        raise ValueError("Cell Ranger count CSV has inconsistent column counts")
    if len(set(header)) != len(header) or any(not key.strip() for key in header):
        raise ValueError("Cell Ranger count CSV has duplicate or blank metric names")
    data = {}
    for key, value in zip(header, values):
        value = value.strip()
        if value.lower() in {"", "na", "n/a", "nan", "none", "null", "-"}:
            continue
        percent = value.endswith("%")
        try:
            number = float(value.removesuffix("%").replace(",", ""))
        except ValueError as error:
            if key in COUNTS | PERCENTAGES:
                raise ValueError(f"Invalid Cell Ranger metric {key!r}: {value!r}") from error
            data[key] = value
            continue
        if percent and key in COUNTS:
            raise ValueError(f"Cell Ranger count metric {key!r} cannot be a percentage")
        if not math.isfinite(number) or number < 0 or (key in PERCENTAGES and number > 100):
            raise ValueError(f"Out-of-range Cell Ranger metric {key!r}: {value!r}")
        data[key] = int(number) if number.is_integer() else number
    if "Estimated Number of Cells" not in data:
        raise ValueError("Cell Ranger count CSV is missing Estimated Number of Cells")
    return data


def sample_name(filename: str, root: str) -> str:
    if filename == "metrics_summary.csv":
        name = Path(root).name
        # Native count layout is <sample>/outs/metrics_summary.csv.
        if name == "outs":
            name = Path(root).parent.name
    else:
        name = filename.removesuffix(".metrics_summary.csv")
    if not name.strip():
        raise ValueError("Cell Ranger count CSV requires a nonblank sample name")
    return name


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Cellranger Count",
            anchor="cellranger_count",
            href="https://www.10xgenomics.com",
            info="10X genomics workflows for single cell sequencing analysis.",
        )
        self.cellranger_count_data = {}
        self._sources = {}
        for f in self.find_log_files("cellranger_count", filehandles=True):
            self.parse_cellranger_count_log(f)
        self.cellranger_count_data = self.ignore_samples(self.cellranger_count_data)
        if not self.cellranger_count_data:
            raise ModuleNoSamplesFound
        log.info("Found %s Cell Ranger count summaries", len(self.cellranger_count_data))
        self.write_data_file(self.cellranger_count_data, "multiqc_cellranger_count")
        self.cellranger_count_set_table_headers()
        self.add_section(
            name="Summary",
            anchor="cellranger_count-summary",
            description="QC metrics from Cell Ranger single cell gene quantification.",
            plot=table.plot(
                self.cellranger_count_data,
                self.cellranger_count_qc_headers,
                {"id": "cellranger_count-summary-table", "title": "Cell Ranger count: Summary"},
            ),
        )
        mapped = self.cellranger_count_mapped_chart()
        if mapped is not None:
            self.add_section(
                name="Confidently Mapped Reads",
                anchor="cellranger_count_confidently_mapped",
                description="Distribution of reads confidently mapped to genomic regions in %. "
                "Samples missing any of the three region metrics are excluded from this plot.",
                plot=mapped,
            )

    def parse_cellranger_count_log(self, f):
        name = sample_name(f["fn"], f["root"])
        source = str(Path(f["root"]) / f["fn"])
        if name in self.cellranger_count_data:
            raise ValueError(f"Duplicate Cell Ranger count sample {name!r}: {self._sources[name]} and {source}")
        try:
            self.cellranger_count_data[name] = parse_summary(f["f"])
        except (ValueError, csv.Error) as error:
            raise ValueError(f"Cannot parse {source}: {error}") from error
        self._sources[name] = source
        self.add_data_source(f, name)

    def cellranger_count_mapped_chart(self):
        data = {
            sample: {key: values[key] for key in REGIONS}
            for sample, values in self.cellranger_count_data.items()
            if all(key in values for key in REGIONS)
        }
        omitted = set(self.cellranger_count_data) - set(data)
        if omitted:
            log.warning(
                "Cell Ranger count region plot excludes samples with missing region metrics: %s",
                ", ".join(sorted(omitted)),
            )
        if not data:
            return None
        return bargraph.plot(
            data,
            REGIONS,
            {
                "id": "cellranger_count_confidently_mapped_plot",
                "title": "Cellranger count: Confidently Mapped Reads",
                "ylab": "% Reads",
                "cpswitch": False,
                "hide_zero_cats": False,
            },
        )

    def cellranger_count_set_table_headers(self):
        count_labels = {
            "Number of Reads": "Reads",
            "Estimated Number of Cells": "# Cells",
            "Mean Reads per Cell": "Mean Reads/Cell",
            "Median UMI Counts per Cell": "Median UMI/Cells",
            "Median Genes per Cell": "Median Genes/Cell",
            "Total Genes Detected": "Detected Genes",
        }
        headers = {
            key: {
                "title": f"{config.read_count_prefix} {title}".strip(),
                "description": key,
                "min": 0,
                "modify": lambda value: value * config.read_count_multiplier,
                "scale": "GnBu",
                "shared_key": "read_count",
            }
            for key, title in count_labels.items()
        }
        percentage_labels = {
            "Valid Barcodes": "% Valid BC",
            "Sequencing Saturation": "Saturation",
            "Q30 Bases in Barcode": "BC % > Q30",
            "Q30 Bases in RNA Read": "RNA % > Q30",
            "Reads Mapped to Genome": "% Mapped to Genome",
            "Reads Mapped Confidently to Genome": "% Mapped Confidently to Genome",
            "Reads Mapped Confidently to Transcriptome": "% Mapped to Transcriptome",
        }
        headers.update(
            {
                key: {
                    "title": title,
                    "description": key,
                    "min": 0,
                    "max": 100,
                    "format": "{:,.1f}",
                    "suffix": "%",
                    "scale": "GnBu",
                }
                for key, title in percentage_labels.items()
            }
        )
        self.cellranger_count_qc_headers = headers
        keys = (
            "Mean Reads per Cell",
            "Estimated Number of Cells",
            "Sequencing Saturation",
            "Reads Mapped Confidently to Genome",
            "Reads Mapped Confidently to Transcriptome",
            "Total Genes Detected",
        )
        self.general_stats_addcols(self.cellranger_count_data, {key: headers[key] for key in keys})
