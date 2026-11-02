"""Test subprocess guard: never connect to SMTP or any network endpoint."""

import sys


def _guard(event, args):
    if event in {"socket.connect", "smtplib.connect"}:
        raise RuntimeError("Network disabled in GCF MultiQC local tests")


sys.addaudithook(_guard)

# Tests supply explicit public/synthetic configuration, never ambient user config.
from multiqc import config  # noqa: E402

config.find_user_files = lambda: None
