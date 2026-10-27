#!/usr/bin/env python3
"""validate_masked_input.py — guard the --masked_genome/--masked_rmout injection path.

When MASK_GENOME derives the masked genome from the same genome file, coordinate
agreement is structural. When a caller injects a pre-masked genome and its
RepeatMasker .out instead, nothing enforces that any more, and every failure
mode is silent: a seqid/assembly mismatch still produces plausible gene counts
and BUSCO scores (AUGUSTUS and HISAT2/StringTie just predict on different
coordinate systems and BUILD_UNION merges them anyway); a seqid mismatch in the
.out reads as "no transposable-element-derived genes in this genome"; a
hard-masked (N-replaced) file just lowers completeness a little. All three look
like normal pipeline output, not an error.

This script is the pre-execution check: exit 1 with the offending seqid(s)
named rather than let any of the three pass silently.

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


def rmout_seqids(path):
    """Seqids named in a RepeatMasker .out (column 5 of each data row). A data
    row is recognised by its first token being the SW score (an integer); this
    skips the fixed 3-line header without depending on its exact wording. An
    empty set is valid (e.g. a genuinely repeat-poor genome's empty .out)."""
    seqids = set()
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                seqids.add(parts[4])
    return seqids


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

    rm_ids = rmout_seqids(rmout_path)
    rm_unknown = sorted(rm_ids - set(genome))
    if rm_unknown:
        errors.append(
            f"{len(rm_unknown)} seqid(s) in --masked_rmout do not appear in --genome at all: "
            + ", ".join(rm_unknown[:10])
            + (" ..." if len(rm_unknown) > 10 else "")
        )

    total_len = sum(v[0] for v in masked.values())
    total_lower = sum(v[1] for v in masked.values())
    total_n = sum(v[2] for v in masked.values())
    lower_frac = total_lower / total_len if total_len else 0.0
    n_frac = total_n / total_len if total_len else 0.0
    if n_frac > 0.01 and lower_frac < 0.001:
        errors.append(
            f"--masked_genome is {n_frac:.1%} N but only {lower_frac:.3%} lowercase: this looks "
            "hard-masked (repeats replaced with N), not soft-masked (repeats lower-cased). "
            "Gene-Miner expects soft-masking; AUGUSTUS and the TE filter both need the repeat "
            "sequence to still be present, only lower-cased."
        )

    if errors:
        print("VALIDATION FAILED for --masked_genome/--masked_rmout:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    print(
        f"masked-genome validation OK: {len(genome)} seqid(s) match between --genome and "
        f"--masked_genome, {len(rm_ids)} seqid(s) referenced in --masked_rmout, "
        f"{lower_frac:.1%} of --masked_genome is lower-case."
    )


if __name__ == "__main__":
    main()
