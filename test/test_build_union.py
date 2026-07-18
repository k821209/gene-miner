"""Smoke tests for bin/build_union.py — the union step that is the heart of the
pipeline (usable AUGUSTUS  UNION  RNA-only  UNION  GeneMark-only, isoform-aware).

Fixtures (all on chr1, + strand):
  AUGUSTUS  g1  score 1.0, 901 bp CDS (300 aa)  -> usable (kept)
            g2  score 0.5, 451 bp CDS           -> fails score>=0.8  (dropped)
            g3  score 1.0, 101 bp CDS (33 aa)    -> fails aa>=100     (dropped)
  RNA-seq   t1  overlaps g1 (extra isoform)   |  t2  no overlap -> RNA-only locus
  GeneMark  gm1 overlaps g1 (extra isoform)   |  gm2 no overlap -> GeneMark-only
            gm3 status "Incomplete"           -> skipped entirely
"""
import os
import tempfile
import unittest

from _helpers import run_script, write, gene_ids

AUG = """chr1\tAUGUSTUS\tgene\t1000\t1900\t1.0\t+\t.\tID=g1
chr1\tAUGUSTUS\tmRNA\t1000\t1900\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tAUGUSTUS\tCDS\t1000\t1900\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
chr1\tAUGUSTUS\tgene\t5000\t5450\t0.5\t+\t.\tID=g2
chr1\tAUGUSTUS\tmRNA\t5000\t5450\t.\t+\t.\tID=g2.t1;Parent=g2
chr1\tAUGUSTUS\tCDS\t5000\t5450\t.\t+\t0\tID=g2.t1.cds;Parent=g2.t1
chr1\tAUGUSTUS\tgene\t7000\t7100\t1.0\t+\t.\tID=g3
chr1\tAUGUSTUS\tmRNA\t7000\t7100\t.\t+\t.\tID=g3.t1;Parent=g3
chr1\tAUGUSTUS\tCDS\t7000\t7100\t.\t+\t0\tID=g3.t1.cds;Parent=g3.t1
"""

TD = """chr1\ttransdecoder\tgene\t1200\t1800\t.\t+\t.\tID=t1
chr1\ttransdecoder\tmRNA\t1200\t1800\t.\t+\t.\tID=t1.m1;Parent=t1
chr1\ttransdecoder\tCDS\t1200\t1500\t.\t+\t0\tID=cds.t1.m1;Parent=t1.m1
chr1\ttransdecoder\tCDS\t1600\t1800\t.\t+\t0\tID=cds.t1.m1;Parent=t1.m1
chr1\ttransdecoder\tgene\t9000\t9600\t.\t+\t.\tID=t2
chr1\ttransdecoder\tmRNA\t9000\t9600\t.\t+\t.\tID=t2.m1;Parent=t2
chr1\ttransdecoder\tCDS\t9000\t9600\t.\t+\t0\tID=cds.t2.m1;Parent=t2.m1
"""

GM = 'chr1\tGeneMark\tCDS\t1300\t1700\t.\t+\t0\tgene_id "gm1"; transcript_id "gm1.t1"; status "complete";\n' \
     'chr1\tGeneMark\tCDS\t12000\t12600\t.\t+\t0\tgene_id "gm2"; transcript_id "gm2.t1"; status "complete";\n' \
     'chr1\tGeneMark\tCDS\t15000\t15600\t.\t+\t0\tgene_id "gm3"; transcript_id "gm3.t1"; status "Incomplete";\n'


def _setup(tmp, with_gm=False):
    write(os.path.join(tmp, "augustus_scaffold.gff3"), AUG)
    write(os.path.join(tmp, "annot", "genome.transdecoder.gff3"), TD)
    args = ["--prefix", "TEST"]
    if with_gm:
        gtf = os.path.join(tmp, "genemark.gtf")
        write(gtf, GM)
        args += ["--genemark", gtf]
    return run_script("build_union.py", args, cwd=tmp)


def _by_prefix(ids):
    return {p: [g for g in ids if "_%s" % p in g] for p in ("A", "R", "E")}


class TestBuildUnion(unittest.TestCase):
    def test_two_stream_union(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = _setup(tmp, with_gm=False)
            self.assertEqual(r.returncode, 0, r.stderr)
            out = os.path.join(tmp, "union.gff3")
            self.assertTrue(os.path.exists(out))
            genes = gene_ids(out)
            groups = _by_prefix(genes)
            # g1 kept as the single AUGUSTUS locus; g2/g3 dropped by QC gates
            self.assertEqual(len(groups["A"]), 1, genes)
            # t2 becomes its own RNA-only locus
            self.assertEqual(len(groups["R"]), 1, genes)
            # no GeneMark stream in this run
            self.assertEqual(len(groups["E"]), 0, genes)
            # the AUGUSTUS locus was augmented with t1's isoform (>= 2 transcripts)
            mrnas = [m for m in gene_ids(out, "mRNA") if groups["A"][0] in m]
            self.assertGreaterEqual(len(mrnas), 2, "AUGUSTUS locus not isoform-augmented")

    def test_three_stream_union_with_genemark(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = _setup(tmp, with_gm=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            out = os.path.join(tmp, "union.gff3")
            groups = _by_prefix(gene_ids(out))
            self.assertEqual(len(groups["A"]), 1)
            self.assertEqual(len(groups["R"]), 1)
            # gm2 -> one GeneMark-only locus; gm1 folded into g1; gm3 skipped (incomplete)
            self.assertEqual(len(groups["E"]), 1, "expected exactly one GeneMark-only locus")
            # g1 now carries its own model + t1 + gm1 -> 3 distinct isoforms
            mrnas = [m for m in gene_ids(out, "mRNA") if groups["A"][0] in m]
            self.assertGreaterEqual(len(mrnas), 3, "GeneMark isoform not folded onto AUGUSTUS locus")

    def test_incomplete_genemark_is_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup(tmp, with_gm=True)
            # gm3 (status Incomplete) must not appear as a locus: only one _E gene (gm2)
            groups = _by_prefix(gene_ids(os.path.join(tmp, "union.gff3")))
            self.assertEqual(len(groups["E"]), 1)


if __name__ == "__main__":
    unittest.main()
