"""Smoke test for bin/lift_agp.py — the contig->scaffold coordinate lift used by
the shell driver (run_gene_mining.sh) for fragmented/RagTag-scaffolded
assemblies. The error-prone case is a reverse-oriented contig, where feature
coordinates are mirrored within the contig and the strand flips (README Lesson 3).

AGP: scaffold1 carries ctgA (1-1000, +) then a gap then ctgB (1101-2100, -).
  ctgA feature 100-200 (+)  -> scaffold1 100-200 (+)      [forward: unchanged here]
  ctgB feature 100-200 (+)  -> scaffold1 1901-2001 (-)    [reverse: mirrored + strand flip]
  ctgC feature (not in AGP) -> dropped
"""
import os
import tempfile
import unittest

from _helpers import run_script, write

AGP = "\n".join([
    "scaffold1\t1\t1000\t1\tW\tctgA\t1\t1000\t+",
    "scaffold1\t1001\t1100\t2\tN\t100\tscaffold\tyes",   # gap line (component_type N) -> ignored
    "scaffold1\t1101\t2100\t3\tW\tctgB\t1\t1000\t-",
]) + "\n"

GFF = "\n".join([
    "##gff-version 3",
    "ctgA\tsrc\tgene\t100\t200\t.\t+\t.\tID=a1",
    "ctgB\tsrc\tgene\t100\t200\t.\t+\t.\tID=b1",
    "ctgC\tsrc\tgene\t50\t100\t.\t+\t.\tID=c1",
]) + "\n"


def _features(gff_path):
    """Return {ID: (seqid, start, end, strand)} for gene lines."""
    import re
    out = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            out[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
    return out


class TestLiftAgp(unittest.TestCase):
    def test_forward_reverse_and_missing_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            agp = os.path.join(tmp, "scaf.agp"); write(agp, AGP)
            gin = os.path.join(tmp, "in.gff3"); write(gin, GFF)
            out = os.path.join(tmp, "lifted.gff3")
            r = run_script("lift_agp.py", [agp, gin, out], cwd=tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("lifted features: 2", r.stdout)
            self.assertIn("skipped: 1", r.stdout)

            feats = _features(out)
            # forward contig: coordinates unchanged (scaffold start = 1), strand kept
            self.assertEqual(feats["a1"], ("scaffold1", 100, 200, "+"))
            # reverse contig: mirrored within a 1000 bp contig placed at 1101, strand flipped
            #   na = 1101 + (1000 - 200) = 1901 ; nb = 1101 + (1000 - 100) = 2001
            self.assertEqual(feats["b1"], ("scaffold1", 1901, 2001, "-"))
            # contig absent from the AGP is dropped
            self.assertNotIn("c1", feats)

    def test_header_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            agp = os.path.join(tmp, "scaf.agp"); write(agp, AGP)
            gin = os.path.join(tmp, "in.gff3"); write(gin, GFF)
            out = os.path.join(tmp, "lifted.gff3")
            run_script("lift_agp.py", [agp, gin, out], cwd=tmp)
            with open(out) as fh:
                self.assertTrue(fh.readline().startswith("##gff-version 3"))


if __name__ == "__main__":
    unittest.main()
