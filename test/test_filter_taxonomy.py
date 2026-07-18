"""Smoke tests for bin/filter_taxonomy.py — the taxonomic QC filter that drops
residual bacterial gene models using eggNOG-mapper orthogroup taxonomy.

Fixtures (emapper.annotations):
  TEST_A000001  OGs rooted only in Bacteria (taxid 2 / 1224)  -> dropped
  TEST_A000002  OGs rooted in Eukaryota / Metazoa (2759/33208) -> kept
  TEST_A000003  no orthogroup ("-")                            -> kept (only 'bac' is dropped)
"""
import os
import tempfile
import unittest

from _helpers import run_script, write, gene_ids, fasta_headers

ANN = "#query\teggNOG_OGs\n" \
      "TEST_A000001\tCOG0001@2|Bacteria,COG0001@1224|Proteobacteria\n" \
      "TEST_A000002\tKOG0001@2759|Eukaryota,KOG0002@33208|Metazoa\n" \
      "TEST_A000003\t-\n"

GFF = """chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=TEST_A000001
chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=TEST_A000001.t1;Parent=TEST_A000001
chr1\tsrc\tCDS\t100\t400\t.\t+\t0\tID=TEST_A000001.t1.cds;Parent=TEST_A000001.t1
chr1\tsrc\tgene\t1000\t1300\t.\t+\t.\tID=TEST_A000002
chr1\tsrc\tmRNA\t1000\t1300\t.\t+\t.\tID=TEST_A000002.t1;Parent=TEST_A000002
chr1\tsrc\tCDS\t1000\t1300\t.\t+\t0\tID=TEST_A000002.t1.cds;Parent=TEST_A000002.t1
chr1\tsrc\tgene\t2000\t2300\t.\t+\t.\tID=TEST_A000003
chr1\tsrc\tmRNA\t2000\t2300\t.\t+\t.\tID=TEST_A000003.t1;Parent=TEST_A000003
chr1\tsrc\tCDS\t2000\t2300\t.\t+\t0\tID=TEST_A000003.t1.cds;Parent=TEST_A000003.t1
"""

PEP = ">TEST_A000001.t1\nMKAILV\n>TEST_A000002.t1\nMQRSTV\n>TEST_A000003.t1\nMHGEDC\n"


class TestFilterTaxonomy(unittest.TestCase):
    def test_bacterial_gene_dropped_eukaryotic_and_unassigned_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = os.path.join(tmp, "emapper.annotations"); write(ann, ANN)
            gin = os.path.join(tmp, "in.gff3"); write(gin, GFF)
            pin = os.path.join(tmp, "in.pep.fa"); write(pin, PEP)
            gout = os.path.join(tmp, "out.gff3")
            pout = os.path.join(tmp, "out.pep.fa")
            r = run_script("filter_taxonomy.py", [ann, gin, pin, gout, pout], cwd=tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("bacterial genes dropped: 1", r.stdout)

            kept = gene_ids(gout)
            self.assertNotIn("TEST_A000001", kept, "bacterial gene not dropped")
            self.assertIn("TEST_A000002", kept, "eukaryotic gene wrongly dropped")
            self.assertIn("TEST_A000003", kept, "unassigned gene wrongly dropped")

            headers = " ".join(fasta_headers(pout))
            self.assertNotIn("TEST_A000001", headers)
            self.assertEqual(len(fasta_headers(pout)), 2)


if __name__ == "__main__":
    unittest.main()
