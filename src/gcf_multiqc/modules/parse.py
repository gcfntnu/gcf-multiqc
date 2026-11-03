"""Parse Biosciences reports, adapted from gcfntnu/MultiQC (GPL-3.0).

The legacy GCF module reports the ``all-sample`` summary for each input library,
not every biological-sample row. Both column- and row-oriented SplitPipe CSVs
are supported without inferring sample identifiers as numbers or missing data.
"""

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

COUNTS = {"number_of_reads", "number_of_cells", "mean_reads_per_cell"}
FRACTIONS = {
    "valid_barcode_fraction",
    "sequencing_saturation",
    "bc1_Q30",
    "bc2_Q30",
    "bc3_Q30",
    "cDNA_Q30",
    "transcriptome_map_fraction",
}
SUFFIXES = ("agg_sample_summary.csv", "agg_samp_ana_summary.csv")


def sample_name(filename: str, root: str) -> str:
    """Use an explicit collector prefix, otherwise the library directory name."""
    for suffix in SUFFIXES:
        if filename == suffix:
            name = Path(root).name
            break
        if filename.endswith("." + suffix):
            name = filename[: -len(suffix) - 1]
            break
    else:
        raise ValueError(f"Unrecognized Parse summary filename: {filename!r}")
    if not name or not name.strip():
        raise ValueError("Parse summary requires a nonblank library name")
    return name


def _value(key: str, value: str):
    value = value.strip()
    if value.lower() in {"", "na", "n/a", "nan", "none", "null", "-"}:
        return None
    if key not in COUNTS | FRACTIONS:
        # Preserve additional producer fields, including textual annotations.
        try:
            number = float(value.replace(",", ""))
        except ValueError:
            return value
        return number if math.isfinite(number) else value
    percent = value.endswith("%")
    if percent and key not in FRACTIONS:
        raise ValueError(f"Count {key!r} cannot be a percentage")
    try:
        number = float(value.removesuffix("%").replace(",", ""))
    except ValueError as error:
        raise ValueError(f"Invalid Parse metric {key!r}: {value!r}") from error
    if percent:
        number /= 100
    if not math.isfinite(number) or number < 0 or (key in FRACTIONS and number > 1):
        raise ValueError(f"Out-of-range Parse metric {key!r}: {value!r}")
    return int(number) if number.is_integer() else number


def parse_summary(handle: TextIO) -> dict:
    rows = [row for row in csv.reader(handle) if row and any(cell.strip() for cell in row)]
    if len(rows) < 2:
        raise ValueError("Parse summary requires a header and data")
    header = rows[0]
    header[0] = header[0].removeprefix("\ufeff")
    if len(header) < 2 or any(len(row) != len(header) for row in rows[1:]):
        raise ValueError("Parse summary has inconsistent CSV column counts")
    if len(set(header[1:])) != len(header[1:]) or any(not key.strip() for key in header[1:]):
        raise ValueError("Parse summary has duplicate or blank column names")
    row_names = [row[0] for row in rows[1:]]
    if len(set(row_names)) != len(row_names) or any(not name.strip() for name in row_names):
        raise ValueError("Parse summary has duplicate or blank row names")
    columns = "all-sample" in header[1:]
    records = "all-sample" in row_names
    if columns == records:
        raise ValueError("Parse summary must contain exactly one all-sample row or column")
    if columns:
        index = header.index("all-sample")
        values = {row[0]: row[index] for row in rows[1:]}
    else:
        values = dict(zip(header[1:], rows[1 + row_names.index("all-sample")][1:]))
    parsed = {key: _value(key, value) for key, value in values.items()}
    parsed = {key: value for key, value in parsed.items() if value is not None}
    if "number_of_reads" not in parsed:
        raise ValueError("Parse all-sample summary is missing number_of_reads")
    return parsed


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Parse Biosciences",
            anchor="parse",
            href="https://www.parsebiosciences.com",
            info="Parse Biosciences workflow for single cell sequencing analysis.",
        )
        self.parse_data = {}
        self._sources = {}
        for f in self.find_log_files("parse", filehandles=True):
            self.parse_log(f)
        self.parse_data = self.ignore_samples(self.parse_data)
        if not self.parse_data:
            raise ModuleNoSamplesFound
        log.info("Found %s Parse summaries", len(self.parse_data))
        self.write_data_file(self.parse_data, "multiqc_parse")
        self.parse_set_table_headers()
        self.add_section(
            name="Summary",
            anchor="parse-summary",
            description="QC metrics from Parse Biosciences single cell workflow.",
            plot=table.plot(
                self.parse_data,
                self.parse_qc_headers,
                {"id": "parse-summary-table", "title": "Parse Biosciences: Summary"},
            ),
        )

    def parse_log(self, f):
        name = sample_name(f["fn"], f["root"])
        source = str(Path(f["root"]) / f["fn"])
        if name in self.parse_data:
            raise ValueError(f"Duplicate Parse library {name!r}: {self._sources[name]} and {source}")
        try:
            self.parse_data[name] = parse_summary(f["f"])
        except (ValueError, csv.Error) as error:
            raise ValueError(f"Cannot parse {source}: {error}") from error
        self._sources[name] = source
        self.add_data_source(f, name)

    def parse_set_table_headers(self):
        labels = {
            "number_of_reads": ("Reads", "Total reads"),
            "number_of_cells": ("# Cells", "Estimated number of cells"),
            "mean_reads_per_cell": ("Mean Reads/Cell", "Mean reads per cell reported by SplitPipe"),
        }
        headers = {
            key: {
                "title": f"{config.read_count_prefix} {title}".strip(),
                "description": description,
                "min": 0,
                "scale": "GnBu",
                "modify": lambda value: value * config.read_count_multiplier,
                "shared_key": "read_count",
            }
            for key, (title, description) in labels.items()
        }
        labels = {
            "valid_barcode_fraction": "% Valid BC",
            "sequencing_saturation": "Saturation",
            "bc1_Q30": "BC1 % > Q30",
            "bc2_Q30": "BC2 % > Q30",
            "bc3_Q30": "BC3 % > Q30",
            "cDNA_Q30": "RNA % > Q30",
            "transcriptome_map_fraction": "% Mapped to Transcriptome",
        }
        headers.update(
            {
                key: {
                    "title": title,
                    "description": title,
                    "min": 0,
                    "max": 100,
                    "modify": lambda value: value * 100,
                    "format": "{:,.1f}",
                    "suffix": "%",
                    "scale": "GnBu",
                }
                for key, title in labels.items()
            }
        )
        self.parse_qc_headers = headers
        keys = (
            "number_of_reads",
            "mean_reads_per_cell",
            "number_of_cells",
            "sequencing_saturation",
            "transcriptome_map_fraction",
        )
        self.general_stats_addcols(self.parse_data, {key: headers[key] for key in keys})
