"""Shared helpers for the gene-miner smoke tests.

The tests exercise the pure-Python catalogue-building logic on tiny synthetic
fixtures. They shell out to the real scripts in ``bin/`` exactly as the pipeline
does, so no heavy dependency (AUGUSTUS, GeneMark-ETP, eggNOG, RepeatMasker,
BUSCO) is needed — only python3. Each test runs in well under a second.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")


def run_script(script, args, cwd=None):
    """Run bin/<script> with the current interpreter; return CompletedProcess."""
    cmd = [sys.executable, os.path.join(BIN, script)] + [str(a) for a in args]
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True) if os.path.dirname(path) else None
    with open(path, "w") as fh:
        fh.write(text)


def gene_ids(gff_path, feature="gene"):
    """IDs of every feature of the given type in a GFF3."""
    import re
    out = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 2 and c[2] == feature:
                m = re.search(r"ID=([^;]+)", c[8])
                if m:
                    out.append(m.group(1))
    return out


def fasta_headers(pep_path):
    with open(pep_path) as fh:
        return [l[1:].strip() for l in fh if l.startswith(">")]
