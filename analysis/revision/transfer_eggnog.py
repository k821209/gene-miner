#!/usr/bin/env python3
"""transfer_eggnog.py — carry eggNOG-mapper rows from a previous catalogue to a
rebuilt one by identical protein sequence, so only genuinely new proteins need
to be re-annotated.

Usage: transfer_eggnog.py old.pep.fa old.emapper.annotations new.pep.fa \
                          out.annotations out.todo.faa
"""
import sys


def fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, ''.join(seq)


old_pep, old_ann, new_pep, out_ann, out_todo = sys.argv[1:6]
seq_of = {n: s.rstrip('*') for n, s in fasta(old_pep)}
row_of = {}
header = []
for line in open(old_ann):
    if line.startswith('#'):
        if line.startswith('#query'):
            header = [line]
        continue
    q = line.split('\t', 1)[0]
    s = seq_of.get(q)
    if s:
        row_of[s] = line.split('\t', 1)[1]

hit = miss = 0
with open(out_ann, 'w') as o, open(out_todo, 'w') as t:
    o.writelines(header)
    for n, s in fasta(new_pep):
        s = s.rstrip('*')
        r = row_of.get(s)
        if r:
            o.write(n + '\t' + r)
            hit += 1
        else:
            t.write(f'>{n}\n{s}\n')
            miss += 1
print(f'transferred {hit} rows, {miss} proteins need eggNOG')
