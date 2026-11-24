"""Smoke tests for bin/merge_into_reference.py -- where an ADDED TRANSCRIPT ends
up, and what that is allowed to do to the reference's own coordinates.

gm_compare.py pairs at the LOCUS level: one best reference gene for the whole
locus. Attaching every transcript of that locus to that gene, and widening the
gene to span them, is wrong in two ways that no check catches, because the file
stays perfectly well-formed:

  * a transcript goes under a gene it barely touches, or does not touch at all;
  * the gene grows over the gene next door and then CONTAINS it, so anything
    asking "which gene holds this feature" gets the neighbour's features too.

Measured on a C. reinhardtii CC-4532 v6.1 build (issue #30): of 26,169 added
transcripts, 888 sat under a gene other than the one their CDS overlaps most,
460 had no CDS overlap with their host at all, and 1,526 were read-through
models filed under one of the two genes they span; 2,272 reference gene spans
were widened. A downstream junction analysis moved 10% of its introns to a
different gene.

Reference layout used throughout -- two ADJACENT genes on the SAME strand,
which is the arrangement that breaks:

  REFA   chr1  1000-2000  +
  REFB   chr1  2100-4000  +
"""
import os
import re
import tempfile
import unittest

from _helpers import run_script, write

REF = """##gff-version 3
chr1\tref\tgene\t1000\t2000\t.\t+\t.\tID=REFA
chr1\tref\tmRNA\t1000\t2000\t.\t+\t.\tID=REFA.1;Parent=REFA
chr1\tref\tCDS\t1200\t1900\t.\t+\t0\tID=REFA.1.cds;Parent=REFA.1
chr1\tref\tgene\t2100\t4000\t.\t+\t.\tID=REFB
chr1\tref\tmRNA\t2100\t4000\t.\t+\t.\tID=REFB.1;Parent=REFB
chr1\tref\tCDS\t2200\t3900\t.\t+\t0\tID=REFB.1.cds;Parent=REFB.1
"""

# Every locus below is PAIRED WITH REFA in the pairs file -- that is the
# locus-level pairing this change stops trusting per transcript.
#
#   GM_1  CDS wholly inside REFB            -> belongs under REFB
#   GM_2  CDS in REFA and in REFB           -> read-through, under neither
#   GM_3  CDS in REFA, span reaches into REFB -> REFA must not be widened
#   GM_4  CDS in REFA, span reaches empty sequence -> REFA may be widened
GM = """##gff-version 3
chr1\tgm\tgene\t2200\t2900\t.\t+\t.\tID=GM_1
chr1\tgm\tmRNA\t2200\t2900\t.\t+\t.\tID=GM_1.t1;Parent=GM_1
chr1\tgm\tCDS\t2200\t2900\t.\t+\t0\tID=GM_1.t1.cds;Parent=GM_1.t1
chr1\tgm\tgene\t1500\t2600\t.\t+\t.\tID=GM_2
chr1\tgm\tmRNA\t1500\t2600\t.\t+\t.\tID=GM_2.t1;Parent=GM_2
chr1\tgm\tCDS\t1500\t1900\t.\t+\t0\tID=GM_2.t1.cds;Parent=GM_2.t1
chr1\tgm\tCDS\t2200\t2600\t.\t+\t0\tID=GM_2.t1.cds;Parent=GM_2.t1
chr1\tgm\tgene\t1500\t2500\t.\t+\t.\tID=GM_3
chr1\tgm\tmRNA\t1500\t2500\t.\t+\t.\tID=GM_3.t1;Parent=GM_3
chr1\tgm\tCDS\t1500\t1900\t.\t+\t0\tID=GM_3.t1.cds;Parent=GM_3.t1
chr1\tgm\tgene\t1500\t1990\t.\t+\t.\tID=GM_4
chr1\tgm\tmRNA\t1500\t1990\t.\t+\t.\tID=GM_4.t1;Parent=GM_4
chr1\tgm\tCDS\t1500\t1900\t.\t+\t0\tID=GM_4.t1.cds;Parent=GM_4.t1
"""

PAIRS = """gm_gene\tref_gene\toverlap\tclass
GM_1\tREFA\t0.60\trevised
GM_2\tREFA\t0.60\trevised
GM_3\tREFA\t0.60\trevised
GM_4\tREFA\t0.60\trevised
"""

# GM_4 alone, so a legitimate widening can be tested without the others moving
# REFA's span first.
GM_WIDEN = """##gff-version 3
chr1\tgm\tgene\t500\t1990\t.\t+\t.\tID=GM_4
chr1\tgm\tmRNA\t500\t1990\t.\t+\t.\tID=GM_4.t1;Parent=GM_4
chr1\tgm\tCDS\t1500\t1900\t.\t+\t0\tID=GM_4.t1.cds;Parent=GM_4.t1
"""
PAIRS_WIDEN = """gm_gene\tref_gene\toverlap\tclass
GM_4\tREFA\t0.60\trevised
"""


def rows(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 8:
                out.append(c)
    return out


def by_id(path):
    return {re.search(r"ID=([^;]+)", c[8]).group(1): c
            for c in rows(path) if "ID=" in c[8]}


def parent_of(path, tid):
    for c in rows(path):
        m = re.search(r"ID=([^;]+)", c[8])
        if m and m.group(1) == tid and c[2] in ("mRNA", "transcript"):
            p = re.search(r"Parent=([^;,]+)", c[8])
            return p.group(1) if p else None
    return None


def attr(c, k):
    m = re.search(rf"(?:^|;){k}=([^;]+)", c[8])
    return m.group(1) if m else None


class _Merge(unittest.TestCase):
    GM_TEXT, PAIRS_TEXT = GM, PAIRS

    def _merge(self, tmp):
        write(os.path.join(tmp, "ref.gff3"), REF)
        write(os.path.join(tmp, "gm.gff3"), self.GM_TEXT)
        write(os.path.join(tmp, "pairs.tsv"), self.PAIRS_TEXT)
        out = os.path.join(tmp, "merged.gff3")
        r = run_script("merge_into_reference.py", [
            "--ref", os.path.join(tmp, "ref.gff3"),
            "--gm", os.path.join(tmp, "gm.gff3"),
            "--pairs", os.path.join(tmp, "pairs.tsv"),
            "--out", out,
        ], cwd=tmp)
        self.assertEqual(r.returncode, 0, f"merge failed:\n{r.stdout}\n{r.stderr}")
        return out, r


class TestPlacedOnItsOwnCDS(_Merge):
    def test_a_transcript_goes_to_the_gene_its_cds_overlaps(self):
        """GM_1 is PAIRED WITH REFA and its CDS is wholly inside REFB. The
        pairing is the locus's; the transcript's own CDS decides."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            self.assertEqual(parent_of(out, "GM_1.t1"), "REFB")

    def test_the_pairing_it_was_moved_off_is_recorded(self):
        """Or the merge report is the only place that fact exists."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            tx = by_id(out)["GM_1.t1"]
            self.assertEqual(attr(tx, "gm_locus_paired_with"), "REFA")
            # 2200..2900 inclusive, against REFB's CDS 2200..3900
            self.assertEqual(attr(tx, "gm_cds_overlap_bp"), "701")


class TestReadThroughGoesUnderNeither(_Merge):
    def test_it_is_not_a_transcript_of_either_gene(self):
        """GM_2's CDS lies in REFA and in REFB. Filing it under one of them is
        the error; it gets a record of its own."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            self.assertNotIn(parent_of(out, "GM_2.t1"), ("REFA", "REFB"))

    def test_its_own_gene_says_what_it_is_and_what_it_spans(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            host = parent_of(out, "GM_2.t1")
            g = by_id(out)[host]
            self.assertEqual(attr(g, "gm_class"), "readthrough")
            self.assertEqual(sorted(attr(g, "gm_readthrough_of").split(",")),
                             ["REFA", "REFB"])
            self.assertEqual(attr(g, "gm_source_locus"), "GM_2")
            # GFF3 requires the parent to span its children.
            self.assertEqual((g[3], g[4]), ("1500", "2600"))

    def test_no_reference_gene_is_attributed_this_locus(self):
        """Not just the Parent: nothing under REFA or REFB may say it came from
        GM_2, or a reader counting a gene's proposals counts this one twice."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            for c in rows(out):
                if c[2] in ("mRNA", "transcript") and \
                        (re.search(r"Parent=(REFA|REFB)\b", c[8])):
                    self.assertNotEqual(attr(c, "gm_source_locus"), "GM_2")


class TestWideningStopsAtTheNeighbour(_Merge):
    def test_the_reference_gene_is_not_widened_over_its_neighbour(self):
        """GM_3's CDS puts it under REFA, but its span reaches 2500 -- inside
        REFB. Widening REFA to 2500 would make REFA contain REFB's CDS."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            a = by_id(out)["REFA"]
            self.assertEqual((a[3], a[4]), ("1000", "2000"))
            self.assertIsNone(attr(a, "gm_original_span"))

    def test_the_transcript_that_needed_it_gets_its_own_record(self):
        """It is not dropped and not clipped -- either would be a quiet lie
        about what was predicted."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            host = parent_of(out, "GM_3.t1")
            self.assertNotIn(host, ("REFA", "REFB"))
            g = by_id(out)[host]
            self.assertEqual(attr(g, "gm_class"), "exceeds_span")
            self.assertEqual(attr(g, "gm_intended_host"), "REFA")

    def test_no_gene_swallows_coding_sequence_of_another(self):
        """The general form, and the harm the issue actually reported: after the
        merge, a reference gene must not contain coding bases of a DIFFERENT
        reference gene that it did not already contain. On the unfixed script
        REFA ran to 2,900 and held REFB's first 701 coding bases, so "which gene
        is this intron in" answered REFA for REFB's introns."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            before, after = by_id(os.path.join(tmp, "ref.gff3")), by_id(out)
            cds = {"REFA": (1200, 1900), "REFB": (2200, 3900)}
            for g in ("REFA", "REFB"):
                for other, (cs, ce) in cds.items():
                    if other == g:
                        continue
                    def held(table):
                        s, e = int(table[g][3]), int(table[g][4])
                        return max(0, min(e, ce) - max(s, cs) + 1)
                    self.assertLessEqual(
                        held(after), held(before),
                        f"{g} gained coding bases of {other}")

    def test_the_refusal_is_reported_not_just_acted_on(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, r = self._merge(tmp)
            self.assertIn("widenings refused (neighbour) 1", r.stdout)
            self.assertIn("REFB", r.stdout)


class TestLegitimateWideningStillHappens(_Merge):
    """The neighbour rule must not become "never widen". GM_4 reaches 500, where
    there is no gene at all, so REFA is widened and says so."""
    GM_TEXT, PAIRS_TEXT = GM_WIDEN, PAIRS_WIDEN

    def test_the_gene_is_widened_into_empty_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, r = self._merge(tmp)
            a = by_id(out)["REFA"]
            self.assertEqual((a[3], a[4]), ("500", "2000"))
            self.assertEqual(attr(a, "gm_original_span"), "1000-2000")
            self.assertEqual(attr(a, "gm_widened_for"), "GM_4")
            self.assertEqual(parent_of(out, "GM_4.t1"), "REFA")
            self.assertIn("widenings refused (neighbour) 0", r.stdout)


class TestNothingIsLost(_Merge):
    def test_every_added_transcript_is_still_in_the_file_exactly_once(self):
        """Four placements, four kinds. A transcript must not be dropped by the
        kind it was sorted into."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            ids = [re.search(r"ID=([^;]+)", c[8]).group(1) for c in rows(out)
                   if c[2] in ("mRNA", "transcript") and c[1] == "gene-miner"]
            self.assertEqual(sorted(ids),
                             ["GM_1.t1", "GM_2.t1", "GM_3.t1", "GM_4.t1"])

    def test_the_reference_is_still_whole(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, r = self._merge(tmp)
            self.assertIn("reference genes in            2", r.stdout)
            self.assertIn("reference genes out           2", r.stdout)
            self.assertIn("duplicate IDs                 0", r.stdout)
            self.assertIn("orphan Parent references      0", r.stdout)
            self.assertIsNone(attr(by_id(out)["REFB.1"], "gm_tier"))


class TestTheAuditorItself(_Merge):
    """bin/check_merge_placement.py asks the same questions from the FILES, so
    the merge's own counts are not the only witness. These hold the auditor
    down -- including the case where it must refuse to answer."""

    def _audit(self, tmp, merged, extra=()):
        return run_script("check_merge_placement.py",
                          ["--ref", os.path.join(tmp, "ref.gff3"),
                           "--merged", merged] + list(extra), cwd=tmp)

    def test_it_passes_the_merge_this_branch_produces(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            r = self._audit(tmp, out)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn("NOT under their largest-CDS-overlap gene        0", r.stdout)

    def test_it_fails_a_file_with_a_transcript_under_the_wrong_gene(self):
        """Hand-built: GM_1's CDS is inside REFB, filed under REFA -- which is
        exactly what the locus-level attachment produced."""
        with tempfile.TemporaryDirectory() as tmp:
            self._merge(tmp)           # writes ref.gff3
            bad = os.path.join(tmp, "bad.gff3")
            write(bad, REF + (
                "chr1\tgene-miner\tmRNA\t2200\t2900\t.\t+\t.\tID=X.t1;Parent=REFA\n"
                "chr1\tgene-miner\tCDS\t2200\t2900\t.\t+\t0\tID=X.t1.cds;Parent=X.t1\n"))
            r = self._audit(tmp, bad)
            self.assertEqual(r.returncode, 1)
            self.assertIn("REFB", r.stdout)

    def test_it_refuses_to_call_a_file_with_no_additions_clean(self):
        """The failure every check here has met once: nothing to look at reads
        as nothing wrong. Exit 2, not 0."""
        with tempfile.TemporaryDirectory() as tmp:
            self._merge(tmp)
            plain = os.path.join(tmp, "plain.gff3")
            write(plain, REF)
            r = self._audit(tmp, plain)
            self.assertEqual(r.returncode, 2)
            self.assertIn("COULD NOT ASSESS", r.stdout)

    def test_a_wrong_source_is_also_could_not_assess(self):
        """Not 'clean'. --source naming a column-2 value the file never uses
        finds no additions, which is indistinguishable from a good merge unless
        it says so."""
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._merge(tmp)
            r = self._audit(tmp, out, ["--source", "not-a-source"])
            self.assertEqual(r.returncode, 2)
            self.assertIn("COULD NOT ASSESS", r.stdout)


if __name__ == "__main__":
    unittest.main()
