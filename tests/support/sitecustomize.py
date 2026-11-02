"""Test subprocess guard: never connect to SMTP or any network endpoint."""
import sys


def _guard(event, args):
    if event in {"socket.connect", "smtplib.connect"}:
        raise RuntimeError("Network disabled in GCF MultiQC local tests")


sys.addaudithook(_guard)
