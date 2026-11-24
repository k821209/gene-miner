"""Smoke tests for bin/merge_into_reference.py — specifically, that the evidence
tier reaches every row a consumer will filter on.

The merged file carries two kinds of addition and they are built differently:

  revised  an extra mRNA under an EXISTING reference gene. The transcript is
           what is being proposed, so gm_tier goes on the mRNA row.
  novel    a new gene plus its own children. The locus is what is being
           proposed, so gm_tier goes on the gene row.

That asymmetry is correct as a description of the proposal, and wrong as an
interface: the obvious filter over a GFF3 is one pass over mRNA rows, and under
it the two kinds answer differently. A novel locus graded `low` would keep all
of its transcripts, silently, as a larger count and no error. So the gene's tier
is repeated on its transcripts, and this test is the filter a consumer writes.

Fixtures:
  reference  REF0001          chr1 1000-2000, one transcript
  ours       GM_R000001       chr1 1000-2000  -> revised, attaches to REF0001
             GM_A000002       chr1 5000-6000  -> novel, tier=low, TWO transcripts
             GM_A000003       chr1 8000-9000  -> novel, tier=high, one transcript
"""
import os
import re
import tempfile
import unittest

from _helpers import run_script, write

REF = """##gff-version 3
chr1\tref\tgene\t1000\t2000\t.\t+\t.\tID=REF0001
chr1\tref\tmRNA\t1000\t2000\t.\t+\t.\tID=REF0001.1;Parent=REF0001
chr1\tref\tCDS\t1000\t2000\t.\t+\t0\tID=REF0001.1.cds;Parent=REF0001.1
"""

# Ours. The novel low-tier locus deliberately has TWO transcripts: the count a
# naive filter gets wrong is the transcript count, not the gene count.
GM = """##gff-version 3
chr1\tgm\tgene\t1000\t2100\t.\t+\t.\tID=GM_R000001
chr1\tgm\tmRNA\t1000\t2100\t.\t+\t.\tID=GM_R000001.t1;Parent=GM_R000001
chr1\tgm\tCDS\t1000\t2100\t.\t+\t0\tID=GM_R000001.t1.cds;Parent=GM_R000001.t1
chr1\tgm\tgene\t5000\t6000\t.\t+\t.\tID=GM_A000002
chr1\tgm\tmRNA\t5000\t6000\t.\t+\t.\tID=GM_A000002.t1;Parent=GM_A000002
chr1\tgm\tCDS\t5000\t6000\t.\t+\t0\tID=GM_A000002.t1.cds;Parent=GM_A000002.t1
chr1\tgm\tmRNA\t5000\t5800\t.\t+\t.\tID=GM_A000002.t2;Parent=GM_A000002
chr1\tgm\tCDS\t5000\t5800\t.\t+\t0\tID=GM_A000002.t2.cds;Parent=GM_A000002.t2
chr1\tgm\tgene\t8000\t9000\t.\t+\t.\tID=GM_A000003
chr1\tgm\tmRNA\t8000\t9000\t.\t+\t.\tID=GM_A000003.t1;Parent=GM_A000003
chr1\tgm\tCDS\t8000\t9000\t.\t+\t0\tID=GM_A000003.t1.cds;Parent=GM_A000003.t1
"""

PAIRS = """gm_gene\tref_gene\toverlap\tclass
GM_R000001\tREF0001\t0.95\trevised
GM_A000002\t\t\tnovel
GM_A000003\t\t\tnovel
"""

# What --tiers supplies: the tier on the gene row, as annotate_tiers.py writes it.
TIERS = """##gff-version 3
chr1\tgm\tgene\t1000\t2100\t.\t+\t.\tID=GM_R000001;tier=high
chr1\tgm\tgene\t5000\t6000\t.\t+\t.\tID=GM_A000002;tier=low
chr1\tgm\tgene\t8000\t9000\t.\t+\t.\tID=GM_A000003;tier=high
"""


def mrna_ids_one_pass(path):
    """The filter a consumer actually writes: one pass, mRNA rows, drop low."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == "mRNA" and "gm_tier=low" not in c[8]:
                m = re.search(r"ID=([^;]+)", c[8])
                if m:
                    out.append(m.group(1))
    return out


def mrna_ids_two_pass(path):
    """The correct answer, found the hard way: collect low-tier gene IDs first,
    then drop a transcript whose own tier is low OR whose parent gene's is."""
    low = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == "gene" and "gm_tier=low" in c[8]:
                m = re.search(r"ID=([^;]+)", c[8])
                if m:
                    low.add(m.group(1))
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == "mRNA" and "gm_tier=low" not in c[8]:
                p = re.search(r"Parent=([^;,]+)", c[8])
                if p and p.group(1) in low:
                    continue
                m = re.search(r"ID=([^;]+)", c[8])
                if m:
                    out.append(m.group(1))
    return out


class TestTierReachesTranscripts(unittest.TestCase):
    def _merge(self, tmp):
        for name, text in (("ref.gff3", REF), ("gm.gff3", GM),
                           ("pairs.tsv", PAIRS), ("tiers.gff3", TIERS)):
            write(os.path.join(tmp, name), text)
        out = os.path.join(tmp, "merged.gff3")
        r = run_script("merge_into_reference.py", [
            "--ref", os.path.join(tmp, "ref.gff3"),
            "--gm", os.path.join(tmp, "gm.gff3"),
            "--pairs", os.path.join(tmp, "pairs.tsv"),
            "--tiers", os.path.join(tmp, "tiers.gff3"),
            "--out", out,
        ], cwd=tmp)
        self.assertEqual(r.returncode, 0, f"merge failed:\n{r.stdout}\n{r.stderr}")
        return out

    def test_one_pass_filter_agrees_with_two_pass(self):
        """The regression this guards: they used to disagree, with no error."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._merge(tmp)
            self.assertEqual(sorted(mrna_ids_one_pass(out)),
                             sorted(mrna_ids_two_pass(out)))

    def test_both_transcripts_of_a_low_novel_locus_are_dropped(self):
        """Not just 'the counts match' — the right transcripts are gone. A novel
        locus can carry several, and each one has to be reachable."""
        with tempfile.TemporaryDirectory() as tmp:
            kept = set(mrna_ids_one_pass(self._merge(tmp)))
            self.assertNotIn("GM_A000002.t1", kept)
            self.assertNotIn("GM_A000002.t2", kept)
            self.assertIn("GM_A000003.t1", kept)      # novel, high
            self.assertIn("GM_R000001.t1", kept)      # revised, high
            self.assertIn("REF0001.1", kept)          # reference rows always stay

    def test_reference_rows_never_gain_a_tier(self):
        """Reference transcripts carry no gm_tier, so they survive any tier
        filter. Writing one onto them would let a filter delete the customer's
        own annotation."""
        with tempfile.TemporaryDirectory() as tmp:
            with open(self._merge(tmp)) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    c = line.rstrip("\n").split("\t")
                    if len(c) > 8 and c[1] != "gene-miner":
                        self.assertNotIn("gm_tier", c[8])

    def test_gene_row_tier_is_unchanged(self):
        """The gene row keeps its tier; repeating it on transcripts adds, never
        moves. gm_compare and the tier tables count genes off this row."""
        with tempfile.TemporaryDirectory() as tmp:
            seen = {}
            with open(self._merge(tmp)) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    c = line.rstrip("\n").split("\t")
                    if len(c) > 8 and c[2] == "gene":
                        i = re.search(r"ID=([^;]+)", c[8])
                        t = re.search(r"gm_tier=([^;]+)", c[8])
                        if i and t:
                            seen[i.group(1)] = t.group(1)
            self.assertEqual(seen.get("GM_A000002"), "low")
            self.assertEqual(seen.get("GM_A000003"), "high")


# --- a reference that types genes with an SO subtype of gene ----------------
#
# Phytozome C. reinhardtii CC-4532 v6.1 writes 810 of its 17,693 genes as
# `transposable_element_gene`. gm_compare.py pairs those correctly, because it
# keys a gene through its mRNAs' Parent rather than through the gene row's
# type; merge_into_reference.py used to index only rows typed `gene`, so the
# pairings it could not resolve were reported as a disagreement about gene IDs
# on a reference that was in fact valid.
#
# `region` is here on purpose: a top-level record with no children. It must NOT
# be counted as a gene, or the counts this script prints become wrong in the
# other direction.
REF_SUBTYPE = """##gff-version 3
chr1\tref\tregion\t1\t20000\t.\t+\t.\tID=chr1_region
chr1\tref\tgene\t1000\t2000\t.\t+\t.\tID=REF0001
chr1\tref\tmRNA\t1000\t2000\t.\t+\t.\tID=REF0001.1;Parent=REF0001
chr1\tref\tCDS\t1000\t2000\t.\t+\t0\tID=REF0001.1.cds;Parent=REF0001.1
chr1\tref\ttransposable_element_gene\t5000\t6000\t.\t+\t.\tID=REFTE001
chr1\tref\tmRNA\t5000\t6000\t.\t+\t.\tID=REFTE001.1;Parent=REFTE001
chr1\tref\tCDS\t5000\t6000\t.\t+\t0\tID=REFTE001.1.cds;Parent=REFTE001.1
"""

# Ours: one locus pairing with the plain gene, one with the subtype gene.
GM_SUBTYPE = """##gff-version 3
chr1\tgm\tgene\t1000\t2000\t.\t+\t.\tID=GM_R000001
chr1\tgm\tmRNA\t1000\t2000\t.\t+\t.\tID=GM_R000001.t1;Parent=GM_R000001
chr1\tgm\tCDS\t1000\t1900\t.\t+\t0\tID=GM_R000001.t1.cds;Parent=GM_R000001.t1
chr1\tgm\tgene\t5000\t6000\t.\t+\t.\tID=GM_R000002
chr1\tgm\tmRNA\t5000\t6000\t.\t+\t.\tID=GM_R000002.t1;Parent=GM_R000002
chr1\tgm\tCDS\t5000\t5900\t.\t+\t0\tID=GM_R000002.t1.cds;Parent=GM_R000002.t1
"""

PAIRS_SUBTYPE = """gm_gene\tref_gene\toverlap\tclass
GM_R000001\tREF0001\t0.95\trevised
GM_R000002\tREFTE001\t0.95\trevised
"""


def feature_ids(path, want):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == want:
                m = re.search(r"ID=([^;]+)", c[8])
                if m:
                    out.append(m.group(1))
    return out


class TestSOGeneSubtypes(unittest.TestCase):
    def _merge(self, tmp):
        for name, text in (("ref.gff3", REF_SUBTYPE), ("gm.gff3", GM_SUBTYPE),
                           ("pairs.tsv", PAIRS_SUBTYPE)):
            write(os.path.join(tmp, name), text)
        out = os.path.join(tmp, "merged.gff3")
        r = run_script("merge_into_reference.py", [
            "--ref", os.path.join(tmp, "ref.gff3"),
            "--gm", os.path.join(tmp, "gm.gff3"),
            "--pairs", os.path.join(tmp, "pairs.tsv"),
            "--out", out,
        ], cwd=tmp)
        return out, r

    def test_a_subtype_gene_can_be_paired(self):
        """The regression: this exited 1 with 'pairing(s) name a reference gene
        not found in --ref', blaming the pairs file for a valid reference."""
        with tempfile.TemporaryDirectory() as tmp:
            out, r = self._merge(tmp)
            self.assertEqual(r.returncode, 0,
                             f"merge failed:\n{r.stdout}\n{r.stderr}")
            self.assertIn("reference genes in            2", r.stdout)
            self.assertIn("reference genes out           2", r.stdout)

    def test_every_reference_gene_goes_in_and_comes_out(self):
        """Both kinds, by ID -- a count can match while the wrong row survives."""
        with tempfile.TemporaryDirectory() as tmp:
            out, r = self._merge(tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(feature_ids(out, "gene"), ["REF0001"])
            self.assertEqual(feature_ids(out, "transposable_element_gene"),
                             ["REFTE001"])

    def test_our_transcript_attaches_to_the_subtype_gene(self):
        """Not merely accepted: the added mRNA has to name the subtype gene as
        its Parent, or the merge silently attached nothing."""
        with tempfile.TemporaryDirectory() as tmp:
            out, r = self._merge(tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            added = []
            with open(out) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    c = line.rstrip("\n").split("\t")
                    if len(c) > 8 and c[1] == "gene-miner" and c[2] == "mRNA":
                        p = re.search(r"Parent=([^;,]+)", c[8])
                        added.append((c[0], p.group(1) if p else None))
            self.assertIn(("chr1", "REFTE001"), added)
            self.assertIn(("chr1", "REF0001"), added)

    def test_a_childless_top_level_record_is_not_a_gene(self):
        """`region` spans the contig and parents nothing. Counting it would
        inflate every gene count this script prints."""
        with tempfile.TemporaryDirectory() as tmp:
            out, r = self._merge(tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("no children -> not a gene record", r.stdout)
            self.assertEqual(len(feature_ids(out, "region")), 1)
            self.assertNotIn("duplicate IDs                 1", r.stdout)


if __name__ == "__main__":
    unittest.main()
