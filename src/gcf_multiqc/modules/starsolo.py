"""STARsolo summary adapter, adapted from gcfntnu/MultiQC (GPL-3.0)."""

from __future__ import annotations

import csv
import logging
import math
from pathlib import Path
from typing import TextIO

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)

COUNTS = {
    "Number of Reads",
    "Estimated Number of Cells",
    "Mean Reads per Cell",
    "Median UMI per Cell",
    "Median Genes per Cell",
    "Total Genes Detected",
}
FRACTIONS = {
    "Reads With Valid Barcodes",
    "Sequencing Saturation",
    "Q30 Bases in CB+UMI",
    "Q30 Bases in RNA read",
    "Reads Mapped to Genome: Unique+Multiple",
    "Reads Mapped to Transcriptome: Unique+Multipe Genes",
    "Reads Mapped to Transcriptome: Unique+Multiple Genes",
}


def parse_summary(handle: TextIO) -> dict:
    data = {}
    for number, row in enumerate(csv.reader(handle), 1):
        if not row or not any(value.strip() for value in row):
            continue
        if len(row) != 2 or not row[0].strip():
            raise ValueError(f"STARsolo row {number} must contain a metric and value")
        key, value = row[0].removeprefix("\ufeff"), row[1].strip()
        if key in data:
            raise ValueError(f"Duplicate STARsolo metric {key!r}")
        # Keep missing metrics explicit while detecting duplicate rows.
        data[key] = None
        if value.lower() in {"", "na", "n/a", "nan", "none", "null", "-"}:
            continue
        percent = value.endswith("%")
        try:
            parsed = float(value.removesuffix("%").replace(",", ""))
        except ValueError as error:
            if key in COUNTS | FRACTIONS:
                raise ValueError(f"Invalid STARsolo metric {key!r}: {value!r}") from error
            data[key] = value
            continue
        if percent:
            if key not in FRACTIONS:
                raise ValueError(f"Unexpected percentage for STARsolo metric {key!r}")
            parsed /= 100
        if not math.isfinite(parsed) or parsed < 0 or (key in FRACTIONS and parsed > 1):
            raise ValueError(f"Out-of-range STARsolo metric {key!r}: {value!r}")
        data[key] = int(parsed) if parsed.is_integer() else parsed
    data = {key: value for key, value in data.items() if value is not None}
    if "Number of Reads" not in data:
        raise ValueError("STARsolo summary is missing Number of Reads")
    return data


def sample_name(filename: str, root: str) -> str:
    name = Path(root).name if filename == "Summary.csv" else filename.removesuffix("_Summary.csv")
    if not name.strip():
        raise ValueError("STARsolo summary requires a nonblank sample name")
    return name


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="starsolo",
            anchor="starsolo",
            href="https://github.com/alexdobin/STAR",
            info="STARsolo single cell gene quantification.",
        )
        self.starsolo_data = {}
        self._sources = {}
        for f in self.find_log_files("starsolo", filehandles=True):
            self.parse_starsolo_log(f)
        self.starsolo_data = self.ignore_samples(self.starsolo_data)
        if not self.starsolo_data:
            raise ModuleNoSamplesFound
        log.info("Found %s STARsolo summaries", len(self.starsolo_data))
        self.write_data_file(self.starsolo_data, "multiqc_starsolo")
        self.starsolo_set_table_headers()
        self.add_section(
            name="STAR solo summary",
            anchor="starsolo-summary",
            description="Mean Reads/Cell uses reads in cells mapped to unique genes divided by estimated cells. "
            "Cell Ranger uses total reads divided by estimated cells, so these metrics differ.",
            plot=table.plot(
                self.starsolo_data,
                self.starsolo_qc_headers,
                {"id": "starsolo-summary-table", "title": "STARsolo: Summary"},
            ),
        )

    def parse_starsolo_log(self, f):
        name = sample_name(f["fn"], f["root"])
        source = str(Path(f["root"]) / f["fn"])
        if name in self.starsolo_data:
            raise ValueError(f"Duplicate STARsolo sample {name!r}: {self._sources[name]} and {source}")
        try:
            self.starsolo_data[name] = parse_summary(f["f"])
        except (ValueError, csv.Error) as error:
            raise ValueError(f"Cannot parse {source}: {error}") from error
        self._sources[name] = source
        self.add_data_source(f, name)

    def starsolo_set_table_headers(self):
        count_labels = {
            "Number of Reads": "Reads",
            "Estimated Number of Cells": "# Cells",
            "Mean Reads per Cell": "Mean Reads/Cell",
            "Median UMI per Cell": "Median UMI/Cells",
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
        fraction_labels = {
            "Reads With Valid Barcodes": "% Valid BC",
            "Sequencing Saturation": "Saturation",
            "Q30 Bases in CB+UMI": "BC % > Q30",
            "Q30 Bases in RNA read": "RNA % > Q30",
            "Reads Mapped to Genome: Unique+Multiple": "% Mapped to Genome",
            "Reads Mapped to Transcriptome: Unique+Multipe Genes": "% Mapped to Transcriptome",
            "Reads Mapped to Transcriptome: Unique+Multiple Genes": "% Mapped to Transcriptome",
        }
        headers.update(
            {
                key: {
                    "title": title,
                    "description": key,
                    "min": 0,
                    "max": 100,
                    "modify": lambda value: value * 100,
                    "format": "{:,.1f}",
                    "suffix": "%",
                    "scale": "GnBu",
                }
                for key, title in fraction_labels.items()
            }
        )
        self.starsolo_qc_headers = headers
        keys = ("Mean Reads per Cell", "Estimated Number of Cells", "Sequencing Saturation")
        self.general_stats_addcols(self.starsolo_data, {key: headers[key] for key in keys})
