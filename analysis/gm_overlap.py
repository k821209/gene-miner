#!/usr/bin/env python3
"""gm_overlap.py — the one definition of the locus-match rule used throughout.

A locus footprint is the MERGED set of its coding bases (the union over its
isoforms, each base counted once). Two loci match when the bases they share
reach a fraction of the smaller footprint:

    ofrac(a, b) = ov(a, b) / min(|a|, |b|)

Accepts either a parsed locus dict with a "cds" key or a bare CDS interval list,
so every script here and bin/gm_compare.py apply the same rule.
"""


def _iv(x):
    return x["cds"] if isinstance(x, dict) else x


def merge(iv):
    out = []
    for s, e in sorted((int(s), int(e)) for s, e in iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def flen(m):
    return sum(e - s + 1 for s, e in m)


def shared(a, b):
    i = j = ov = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi >= lo:
            ov += hi - lo + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return ov


def ofrac(a, b):
    ma, mb = merge(_iv(a)), merge(_iv(b))
    return shared(ma, mb) / (min(flen(ma), flen(mb)) or 1)
