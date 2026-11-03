"""Legacy 10x mkfastq JSON adapter from gcfntnu/MultiQC (GPL-3.0).

This is the GCF meaning of ``cellranger`` (qc_summary.json), distinct from the
upstream count/VDJ web-summary parser. The registration hook owns that mapping.
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import TextIO

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)

COUNTS = {"number_reads", "gem_count_estimate"}
FRACTIONS = {
    "barcode_exact_match_ratio",
    "barcode_q30_base_ratio",
    "read1_q30_base_ratio",
    "read2_q30_base_ratio",
    "bc_on_whitelist",
}


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate mkfastq JSON key {key!r}")
        result[key] = value
    return result


def parse_summary(handle: TextIO) -> tuple[str, dict]:
    data = json.load(handle, object_pairs_hook=_unique_object)
    if not isinstance(data, dict):
        raise ValueError("mkfastq summary must be a JSON object")
    version = data.get("10x_software_version")
    samples = data.get("sample_qc")
    if not isinstance(version, str) or not version.strip():
        raise ValueError("mkfastq summary requires 10x_software_version")
    if not isinstance(samples, dict):
        raise ValueError("mkfastq summary requires a sample_qc object")
    result = {}
    for name, sample in samples.items():
        if not name.strip() or not isinstance(sample, dict) or not isinstance(sample.get("all"), dict):
            raise ValueError(f"mkfastq sample {name!r} requires a nonblank name and an all metrics object")
        metrics = dict(sample["all"])
        for key in COUNTS | FRACTIONS | {"mean_barcode_qscore"}:
            if key not in metrics:
                continue
            raw = metrics[key]
            if raw is None or (isinstance(raw, str) and raw.strip().lower() in {"", "na", "n/a", "none", "null", "-"}):
                del metrics[key]
                continue
            if isinstance(raw, bool) or not isinstance(raw, (int, float, str)):
                raise ValueError(f"Invalid mkfastq metric {key!r} for {name!r}: {raw!r}")
            value = str(raw).strip()
            percent = value.endswith("%")
            if percent and key not in FRACTIONS:
                raise ValueError(f"mkfastq count/score {key!r} cannot be a percentage")
            try:
                number = float(value.removesuffix("%").replace(",", ""))
            except ValueError as error:
                raise ValueError(f"Invalid mkfastq metric {key!r} for {name!r}: {raw!r}") from error
            if percent:
                number /= 100
            if not math.isfinite(number) or number < 0 or (key in FRACTIONS and number > 1):
                raise ValueError(f"Out-of-range mkfastq metric {key!r} for {name!r}: {raw!r}")
            metrics[key] = int(number) if number.is_integer() else number
        if "number_reads" not in metrics:
            raise ValueError(f"mkfastq sample {name!r} is missing number_reads")
        result[name] = metrics
    return version, result


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # The new explicit alias can coexist with upstream Cell Ranger. Only the
        # historical selector retains the old anchor and export filename.
        explicit_alias = self.mod_id == "gcf_cellranger_mkfastq"
        anchor = "gcf_cellranger_mkfastq" if explicit_alias else "cellranger"
        super().__init__(
            name="GCF mkfastq" if explicit_alias else "cellranger",
            anchor=anchor,
            href="https://www.10xgenomics.com",
            info="10X genomics mkfastq sequencing and demultiplexing metrics.",
        )
        self.cellranger_data = {}
        self._sources = {}
        self.cellranger_versions = set()
        for f in self.find_log_files(f"{self.mod_id or 'cellranger'}/mkfastq", filehandles=True):
            self.parse_cellranger_log(f)
        self.cellranger_data = self.ignore_samples(self.cellranger_data)
        if not self.cellranger_data:
            raise ModuleNoSamplesFound
        log.info("Found %s mkfastq sample summaries", len(self.cellranger_data))
        self.write_data_file(self.cellranger_data, f"multiqc_{anchor}")
        self.cellranger_set_table_headers()
        self.add_section(
            name="10X genomics mkfastq",
            anchor=f"{anchor}-mkfastq-qc",
            description="QC metrics from 10X genomics mkfastq pipelines "
            "(Cell Ranger, Space Ranger, Cell Ranger ATAC and Long Ranger).",
            plot=table.plot(
                self.cellranger_data,
                self.cellranger_qc_headers,
                {"id": f"{anchor}-mkfastq-qc-table", "title": "10X genomics mkfastq: QC"},
            ),
        )

    def parse_cellranger_log(self, f):
        source = str(Path(f["root"]) / f["fn"])
        try:
            version, samples = parse_summary(f["f"])
        except (ValueError, TypeError) as error:
            raise ValueError(f"Cannot parse {source}: {error}") from error
        for name in samples:
            if name in self.cellranger_data:
                raise ValueError(f"Duplicate mkfastq sample {name!r}: {self._sources[name]} and {source}")
        self.cellranger_versions.add(version)
        self.cellranger_version = version  # Legacy attribute, not used for mixed-run unit labels.
        for name, metrics in samples.items():
            self.cellranger_data[name] = metrics
            self._sources[name] = source
            self.add_data_source(f, name)
            self.add_software_version(version, sample=name, software_name="10x mkfastq")

    def cellranger_set_table_headers(self):
        atac_only = all("atac" in version.lower() for version in self.cellranger_versions)
        title = "PE Reads" if atac_only else "Reads"
        headers = {
            "number_reads": {
                "title": f"{config.read_count_prefix} {title}".strip(),
                "description": "Total reads before filtering, as reported by mkfastq. "
                "Cell Ranger ATAC reports paired-end counts.",
                "min": 0,
                "modify": lambda value: value * config.read_count_multiplier,
                "scale": "GnBu",
                "shared_key": "read_count",
            }
        }
        if any("gem_count_estimate" in values for values in self.cellranger_data.values()):
            headers["gem_count_estimate"] = {
                "title": f"{config.read_count_prefix} GEM estimate".strip(),
                "description": "GEM count estimate",
                "min": 0,
                "modify": lambda value: value * config.read_count_multiplier,
                "scale": "GnBu",
                "shared_key": "read_count",
            }
        labels = {
            "barcode_exact_match_ratio": "Ratio perfect BC",
            "bc_on_whitelist": "BC on whitelist",
            "barcode_q30_base_ratio": "BC % > Q30",
            "read1_q30_base_ratio": "R1 % > Q30",
            "read2_q30_base_ratio": "R2 % > Q30",
        }
        headers.update(
            {
                key: {
                    "title": label,
                    "description": label,
                    "min": 0,
                    "max": 100,
                    "modify": lambda value: value * 100,
                    "format": "{:,.1f}",
                    "suffix": "%",
                    "scale": "GnBu",
                }
                for key, label in labels.items()
            }
        )
        headers["mean_barcode_qscore"] = {
            "title": "BC Qscore",
            "description": "Mean barcode Q score",
            "format": "{:,.1f}",
            "min": 0,
            "scale": "GnBu",
        }
        self.cellranger_qc_headers = headers
        # The legacy mkfastq module intentionally adds no General Statistics columns.
