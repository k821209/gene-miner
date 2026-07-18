"""Smoke tests for bin/te_filter_genes.py — the TE (transposable-element) QC
filter that drops genes whose CDS lies mostly inside interspersed repeats.

Fixtures:
  RepeatMasker .out : chr1 1000-2000 LINE/L1 (interspersed)  -> counts
                      chr1 5000-6000 Simple_repeat            -> ignored (SKIP class)
  genes  TEST_A000001  CDS 1100-1900  fully inside the LINE   -> dropped
         TEST_A000002  CDS 3000-3500  no repeat               -> kept
         TEST_A000003  CDS 5100-5900  inside a Simple_repeat  -> kept (simple repeats
                                                                  are not TE evidence)
"""
import os
import tempfile
import unittest

from _helpers import run_script, write, gene_ids, fasta_headers

# RepeatMasker .out: a 3-line header (skipped: col0 not numeric) then data rows.
RM_OUT = """   SW   perc perc perc  query      position in query           matching repeat
score  div. del. ins.  sequence   begin  end    (left)   repeat      class/family

 1200  5.0  0.0  0.0  chr1        1000   2000   (8000)  + L1elem      LINE/L1        1  1001  (0)  1
  500 10.0  0.0  0.0  chr1        5000   6000   (4000)  + (AT)n       Simple_repeat  1  1001  (0)  2
"""

GFF = """chr1\tsrc\tgene\t1100\t1900\t.\t+\t.\tID=TEST_A000001
chr1\tsrc\tmRNA\t1100\t1900\t.\t+\t.\tID=TEST_A000001.t1;Parent=TEST_A000001
chr1\tsrc\tCDS\t1100\t1900\t.\t+\t0\tID=TEST_A000001.t1.cds;Parent=TEST_A000001.t1
chr1\tsrc\tgene\t3000\t3500\t.\t+\t.\tID=TEST_A000002
chr1\tsrc\tmRNA\t3000\t3500\t.\t+\t.\tID=TEST_A000002.t1;Parent=TEST_A000002
chr1\tsrc\tCDS\t3000\t3500\t.\t+\t0\tID=TEST_A000002.t1.cds;Parent=TEST_A000002.t1
chr1\tsrc\tgene\t5100\t5900\t.\t+\t.\tID=TEST_A000003
chr1\tsrc\tmRNA\t5100\t5900\t.\t+\t.\tID=TEST_A000003.t1;Parent=TEST_A000003
chr1\tsrc\tCDS\t5100\t5900\t.\t+\t0\tID=TEST_A000003.t1.cds;Parent=TEST_A000003.t1
"""

PEP = ">TEST_A000001.t1\nMKAILV\n>TEST_A000002.t1\nMQRSTV\n>TEST_A000003.t1\nMHGEDC\n"


class TestTEFilter(unittest.TestCase):
    def test_te_gene_dropped_others_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            rm = os.path.join(tmp, "rm.out"); write(rm, RM_OUT)
            gin = os.path.join(tmp, "in.gff3"); write(gin, GFF)
            pin = os.path.join(tmp, "in.pep.fa"); write(pin, PEP)
            gout = os.path.join(tmp, "out.gff3")
            pout = os.path.join(tmp, "out.pep.fa")
            r = run_script("te_filter_genes.py", [rm, gin, pin, gout, pout, "0.5"], cwd=tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("TE-derived genes dropped: 1", r.stdout)

            kept = gene_ids(gout)
            self.assertNotIn("TEST_A000001", kept, "TE gene not dropped")
            self.assertIn("TEST_A000002", kept, "clean gene wrongly dropped")
            self.assertIn("TEST_A000003", kept, "simple-repeat gene wrongly treated as TE")

            headers = " ".join(fasta_headers(pout))
            self.assertNotIn("TEST_A000001", headers)
            self.assertIn("TEST_A000002", headers)
            self.assertIn("TEST_A000003", headers)

    def test_threshold_is_respected(self):
        # With a threshold above the observed overlap fraction, nothing is dropped.
        with tempfile.TemporaryDirectory() as tmp:
            rm = os.path.join(tmp, "rm.out"); write(rm, RM_OUT)
            gin = os.path.join(tmp, "in.gff3"); write(gin, GFF)
            pin = os.path.join(tmp, "in.pep.fa"); write(pin, PEP)
            gout = os.path.join(tmp, "out.gff3")
            pout = os.path.join(tmp, "out.pep.fa")
            r = run_script("te_filter_genes.py", [rm, gin, pin, gout, pout, "1.01"], cwd=tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("TE-derived genes dropped: 0", r.stdout)
            self.assertIn("TEST_A000001", gene_ids(gout))


if __name__ == "__main__":
    unittest.main()
