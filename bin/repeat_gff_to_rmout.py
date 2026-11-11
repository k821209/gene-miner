#!/usr/bin/env python3
"""repeat_gff_to_rmout.py -- turn a published repeat-annotation GFF into the two
files Gene-Miner's injection path wants: a soft-masked genome and a
RepeatMasker .out. Driven by --repeat_gff; main.nf then runs its own
VALIDATE_MASKED_INPUT gate over the result.

Why it is in the pipeline instead of left to the caller: published repeat
annotations come as GFF, not as a RepeatMasker .out, and the conversion has
three traps that are each a silent wrong answer rather than an error.

1. The query "(left)" column must be `seqlen - end` computed from the genome
   being masked. RepeatMasker defines it as the number of bases AFTER the match
   on that sequence, and bin/validate_masked_input.py reconstructs the sequence
   length as end + left to catch a .out belonging to a different assembly
   version. Writing left from anything else either defeats that check or trips
   it on a file that is correct.
2. Column 11 -- the repeat class -- decides results. te_filter_genes.py skips
   Simple_repeat / Low_complexity / Satellite on the grounds that real genes
   carry those. Writing every row "Unknown" therefore makes the TE filter treat
   every simple repeat as an interspersed repeat and discard genuine genes whose
   CDS contains a microsatellite, as a lower gene count with no error. On one
   published GFF that was 18% of all intervals (152,144 of 844,228).
3. RepeatMasker writes GFF2 attributes -- Target "Motif:(AT)n" -- not the GFF3
   Target=. A GFF3 parser silently labels every interval alike, which is how (2)
   happens.

Worked example (SGN ships these beside the annotation, on the assembly's own
seqids): ITAG4.0_RepeatModeler_repeats_light.gff, 844,228 intervals, converted
to 63.4% soft-masked and classed Unknown 692,084 / Simple_repeat 126,137 /
Low_complexity 26,007. Using it skips de-novo RepeatModeler, the most expensive
step in the pipeline (8-9 h per multi-Gb genome here).

Usage: repeat_gff_to_rmout.py <genome.fa[.gz]> <repeats.gff[.gz]> <out_prefix>
Writes <out_prefix>.softmasked.fa and <out_prefix>.out
"""
import sys, re, gzip, collections

def op(p, mode="rt"):
    return gzip.open(p, mode, errors="replace") if p.endswith(".gz") else open(p, errors="replace")


# The class column is NOT decoration: te_filter_genes.py reads column 11 and
# SKIPS Simple_repeat / Low_complexity / Satellite, on the grounds that real
# genes carry those. Writing every row as "Unknown" would therefore make the
# TE filter treat ~54,000 simple-repeat intervals as interspersed repeats and
# throw out genuine genes whose CDS happens to contain a microsatellite --
# silently, as a lower gene count.
SIMPLE = re.compile(r"^\(.+\)n$")                       # (AT)n, (TA)n, (T)n ...
LOWCOMP = re.compile(r"-rich$|^poly(purine|pyrimidine)", re.I)   # A-rich, GA-rich ...

def classify(name):
    if SIMPLE.match(name):
        return "Simple_repeat"
    if LOWCOMP.search(name):
        return "Low_complexity"
    return "Unknown"            # rnd-N_family-N: de-novo families, genuinely
                                # unclassified but genuinely interspersed


def repeat_name(attr):
    """RepeatMasker writes GFF2 here: Target "Motif:rnd-4_family-195" 1 6592.
    The GFF3 form (Target=...;) does not appear, so a Target= parser silently
    labels every interval Unspecified -- which is how the class above was lost."""
    m = re.search(r'Target\s+"Motif:([^"]+)"', attr)
    if m:
        return m.group(1)
    m = re.search(r'(?:Target|Name|ID)=([^;\s]+)', attr)
    return m.group(1) if m else "Unspecified"

HEADER = (
    "   SW   perc perc perc  query      position in query           matching       "
    "repeat              position in  repeat\n"
    " score   div. del. ins.  sequence   begin  end    (left)   repeat         "
    "class/family    begin  end (left)   ID\n\n"
)

def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    fa, gff, pref = sys.argv[1:4]

    # pass 1: sequence lengths, needed for (left) and for the mask pass
    order, lens = [], {}
    sid = None
    with op(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                sid = line[1:].split()[0]; order.append(sid); lens[sid] = 0
            else:
                lens[sid] += len(line.strip())
    print(f"FASTA: {len(order)} seqids, {sum(lens.values()):,} bp")

    # pass 2: repeat intervals per seqid
    iv = collections.defaultdict(list)
    skipped = collections.Counter()
    with op(gff) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            s = f[0]
            if s not in lens:
                skipped[s] += 1
                continue
            try:
                b, e = int(f[3]), int(f[4])
            except ValueError:
                continue
            if e > lens[s]:                      # clamp rather than emit an interval
                e = lens[s]                      # the guard would (correctly) reject
            if b < 1:
                b = 1
            if b > e:
                continue
            strand = f[6] if f[6] in "+-" else "+"
            name = repeat_name(f[8])
            try:
                div = float(f[5])
            except (ValueError, IndexError):
                div = 0.0
            iv[s].append((b, e, strand, name, classify(name), div))
    if skipped:
        print(f"  WARNING: {sum(skipped.values()):,} GFF rows on {len(skipped)} seqid(s) "
              f"absent from the FASTA, dropped: {', '.join(list(skipped)[:5])}")
    total = sum(len(v) for v in iv.values())
    cls_n = collections.Counter(r[4] for v in iv.values() for r in v)
    print(f"GFF  : {total:,} intervals on {len(iv)} seqids")
    print("  class (te_filter_genes.py skips Simple_repeat / Low_complexity): "
          + ", ".join(f"{k}={v:,}" for k, v in cls_n.most_common()))

    # .out — one data row per interval, (left) computed from THIS fasta
    with open(f"{pref}.out", "w") as out:
        out.write(HEADER)
        rid = 0
        for s in order:
            for b, e, strand, name, cls, div in sorted(iv.get(s, [])):
                rid += 1
                left = lens[s] - e                       # the invariant the guard checks
                rlen = e - b + 1
                out.write(f"  255   {div:4.1f}  0.0  0.0  {s}  {b}  {e} ({left}) "
                          f"{'+' if strand=='+' else 'C'}  {name}  {cls}  "
                          f"1  {rlen} (0)  {rid}\n")
    print(f"wrote {pref}.out  ({rid:,} data rows)")

    # soft-masked FASTA — lower-case inside intervals, sequence preserved
    masked_bp = 0
    with op(fa) as fh, open(f"{pref}.softmasked.fa", "w") as out:
        sid, buf = None, []
        def flush():
            nonlocal masked_bp
            if sid is None:
                return
            seq = list("".join(buf))
            for b, e, _, _, _, _ in iv.get(sid, []):
                for i in range(b - 1, min(e, len(seq))):
                    if not seq[i].islower():
                        masked_bp += 1
                    seq[i] = seq[i].lower()
            out.write(f">{sid}\n")
            s = "".join(seq)
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + "\n")
        for line in fh:
            if line.startswith(">"):
                flush(); sid = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        flush()
    tot = sum(lens.values())
    print(f"wrote {pref}.softmasked.fa  (soft-masked {masked_bp:,} bp = {masked_bp/tot*100:.2f}%)")

if __name__ == "__main__":
    main()
