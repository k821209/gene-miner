#!/usr/bin/env python3
"""validate_masked_input.py — guard the --masked_genome/--masked_rmout injection path.

When MASK_GENOME derives the masked genome from the same genome file, coordinate
agreement is structural. When a caller injects a pre-masked genome and its
RepeatMasker .out instead, nothing enforces that any more, and every failure
mode is silent: a seqid/assembly mismatch still produces plausible gene counts
and BUSCO scores (AUGUSTUS and HISAT2/StringTie just predict on different
coordinate systems and BUILD_UNION merges them anyway); a seqid mismatch in the
.out reads as "no transposable-element-derived genes in this genome"; an
unmasked or hard-masked file just lowers completeness a little. All of these
look like normal pipeline output, not an error.

Same-named seqid, different assembly version is the sharpest case: the seqid
check alone passes (the name matches), but the .out's own coordinates can
still belong to a different-length chromosome. RepeatMasker's query "(left)"
column is seqlen - end, so end + left reconstructs the chromosome length the
.out was generated against; comparing that to --genome's actual length catches
a version mismatch a name-only check cannot.

This script is the pre-execution check: exit 1 with the offending seqid(s)
named rather than let any of these pass silently.

Usage: validate_masked_input.py genome.fa masked_genome.fa masked_rmout.out
"""
import sys


def fasta_lengths(path):
    """seqid -> (length, lowercase_count, n_count), from a single streaming pass."""
    lengths = {}
    seqid = None
    length = lower = n_count = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seqid is not None:
                    lengths[seqid] = (length, lower, n_count)
                seqid = line[1:].split()[0]
                length = lower = n_count = 0
            else:
                s = line.strip()
                length += len(s)
                lower += sum(1 for c in s if c.islower())
                n_count += sum(1 for c in s if c in "Nn")
        if seqid is not None:
            lengths[seqid] = (length, lower, n_count)
    return lengths


def rmout_rows(path):
    """Yield (seqid, end, left) per RepeatMasker .out data row. A data row is
    recognised by its first token being the SW score (an integer); this skips
    the fixed 3-line header without depending on its exact wording. `left` is
    the query "(left)" column (seqlen - end) with its parentheses stripped;
    rows where it cannot be parsed as an int are skipped rather than crashing,
    since a row this script cannot read cannot be checked either way."""
    for line in open(path):
        parts = line.split()
        if len(parts) >= 8 and parts[0].isdigit():
            seqid = parts[4]
            try:
                end = int(parts[6])
                left = int(parts[7].strip("()"))
            except ValueError:
                continue
            yield seqid, end, left


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: validate_masked_input.py genome.fa masked_genome.fa masked_rmout.out")
    genome_path, masked_path, rmout_path = sys.argv[1:4]

    genome = fasta_lengths(genome_path)
    masked = fasta_lengths(masked_path)

    errors = []

    missing_in_masked = sorted(set(genome) - set(masked))
    extra_in_masked = sorted(set(masked) - set(genome))
    if missing_in_masked:
        errors.append(
            f"{len(missing_in_masked)} seqid(s) in --genome are absent from --masked_genome: "
            + ", ".join(missing_in_masked[:10])
            + (" ..." if len(missing_in_masked) > 10 else "")
        )
    if extra_in_masked:
        errors.append(
            f"{len(extra_in_masked)} seqid(s) in --masked_genome are absent from --genome "
            "(different assembly/version?): "
            + ", ".join(extra_in_masked[:10])
            + (" ..." if len(extra_in_masked) > 10 else "")
        )

    len_mismatches = []
    for sid in sorted(set(genome) & set(masked)):
        g_len = genome[sid][0]
        m_len = masked[sid][0]
        if g_len != m_len:
            len_mismatches.append(f"{sid} ({g_len} bp in --genome vs {m_len} bp in --masked_genome)")
    if len_mismatches:
        errors.append(
            f"{len(len_mismatches)} seqid(s) have a different length between --genome and "
            "--masked_genome (coordinate systems do not match): " + "; ".join(len_mismatches[:10])
            + (" ..." if len(len_mismatches) > 10 else "")
        )

    # Seqid presence in the .out is necessary but not sufficient: a same-named
    # seqid from a different assembly version still passes a name check, so
    # reconstruct the chromosome length each row's coordinates imply
    # (end + left) and compare it to --genome's actual length for that seqid.
    rm_ids = set()
    implied_len = {}       # seqid -> one reconstructed length already confirmed (for internal-consistency check)
    version_mismatches = []
    internal_mismatches = []
    for seqid, end, left in rmout_rows(rmout_path):
        rm_ids.add(seqid)
        this_len = end + left
        if seqid in genome:
            g_len = genome[seqid][0]
            if this_len != g_len and seqid not in [m.split()[0] for m in version_mismatches]:
                version_mismatches.append(
                    f"{seqid} (a repeat at query {end} bp implies a {this_len} bp chromosome via "
                    f"end+left, but --genome's {seqid} is {g_len} bp)"
                )
        else:
            prev = implied_len.get(seqid)
            if prev is None:
                implied_len[seqid] = this_len
            elif prev != this_len and seqid not in [m.split()[0] for m in internal_mismatches]:
                internal_mismatches.append(f"{seqid} (rows imply both {prev} bp and {this_len} bp)")

    rm_unknown = sorted(rm_ids - set(genome))
    if rm_unknown:
        errors.append(
            f"{len(rm_unknown)} seqid(s) in --masked_rmout do not appear in --genome at all: "
            + ", ".join(rm_unknown[:10])
            + (" ..." if len(rm_unknown) > 10 else "")
        )
    if version_mismatches:
        errors.append(
            f"{len(version_mismatches)} seqid(s) in --masked_rmout have the same name as a "
            "--genome seqid but imply a different chromosome length from their own coordinates "
            "(same seqid, different assembly version): " + "; ".join(version_mismatches[:5])
            + (" ..." if len(version_mismatches) > 5 else "")
        )
    if internal_mismatches:
        errors.append(
            "--masked_rmout is internally inconsistent, implying more than one length for the "
            "same seqid not present in --genome: " + "; ".join(internal_mismatches[:5])
        )

    total_len = sum(v[0] for v in masked.values())
    total_lower = sum(v[1] for v in masked.values())
    total_n = sum(v[2] for v in masked.values())
    lower_frac = total_lower / total_len if total_len else 0.0
    n_frac = total_n / total_len if total_len else 0.0
    # A caller asserting --masked_genome IS masked should see real soft-masking;
    # unlike MASK_GENOME's own de-novo path, there is no legitimate reason for an
    # injected "masked" genome to show essentially none (checked independently of
    # n_frac, so an unmasked raw genome fails here too, not just a hard-masked one).
    if lower_frac < 0.001:
        if n_frac > 0.01:
            reason = f"looks hard-masked (repeats replaced with N, {n_frac:.1%} of bases) rather than soft-masked"
        else:
            reason = "looks entirely unmasked (no lower-case, no N)"
        errors.append(
            f"--masked_genome is only {lower_frac:.3%} lower-case: {reason}. Gene-Miner expects "
            "soft-masking; AUGUSTUS and the TE filter both need the repeat sequence to still be "
            "present, only lower-cased."
        )

    if errors:
        print("VALIDATION FAILED for --masked_genome/--masked_rmout:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    print(
        f"masked-genome validation OK: {len(genome)} seqid(s) match between --genome and "
        f"--masked_genome (lengths and RepeatMasker end+left coordinates agree), "
        f"{len(rm_ids)} seqid(s) referenced in --masked_rmout, "
        f"{lower_frac:.1%} of --masked_genome is lower-case."
    )


if __name__ == "__main__":
    main()
