"""Preserve literal custom-content sample IDs while using the upstream renderer."""

from multiqc import config
from multiqc.core.special_case_modules.custom_content import custom_module_classes as upstream_custom_content


def custom_module_classes():
    # GCF's historical custom-content reader did not clean sample-key dictionaries.
    # Scope this compatibility behavior to custom content; standard tool parsers
    # still need filename cleaning to join their measurements to the same sample.
    clean_names = config.fn_clean_sample_names
    config.fn_clean_sample_names = False
    try:
        return upstream_custom_content()
    finally:
        config.fn_clean_sample_names = clean_names
