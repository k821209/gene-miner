#!/usr/bin/env python3
import re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib_venn import venn3, venn3_circles
try:
    from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm
    EQ = DefaultLayoutAlgorithm(fixed_subset_sizes=(1,1,1,1,1,1,1))
except Exception:
    EQ = None
GREEN,BLUE,ORANGE='#4e9a5f','#5b9bd5','#e8983a'
# Counts are read from revised_features.py output (representative-transcript unit),
# so the figure cannot drift from the tables: pass the two files on the command line.
import sys
def load(path, crop):
    txt = open(path).read()
    rep = txt.split("-- representative")[1].split("-- aggregate")[0]
    total = int(re.search(r"revised_loci\t(\d+)", txt).group(1))
    v = {tuple(x == "True" for x in k.split(", ")): int(n)
         for k, n in re.findall(r"\((True|False, (?:True|False), (?:True|False))\):(\d+)", rep.replace("(True", "(True").replace("(False", "(False"))}
    v = {}
    for k, n in re.findall(r"\(((?:True|False), (?:True|False), (?:True|False))\):(\d+)", rep):
        v[tuple(x == "True" for x in k.split(", "))] = int(n)
    T, F = True, False
    sub = (v.get((T,F,F),0), v.get((F,T,F),0), v.get((T,T,F),0), v.get((F,F,T),0),
           v.get((T,F,T),0), v.get((F,T,T),0), v.get((T,T,T),0))
    other = v.get((F,F,F),0)
    iso = sum(n for k,n in v.items() if k[0]); exon = sum(n for k,n in v.items() if k[1]); cds = sum(n for k,n in v.items() if k[2])
    f = lambda x: f"{x:,}"
    return dict(crop=crop, total=f(total), circ=f(total-other), other=f(other), subsets=sub,
                iso=f(iso), exon=f(exon), cds=f(cds))
P = {"a": load(sys.argv[1], "rice"), "b": load(sys.argv[2], "soybean")}
OUT = sys.argv[3] if len(sys.argv) > 3 else "/tmp/fig2_new"
fig,axes=plt.subplots(1,2,figsize=(13,6.8))
for key,ax in zip(('a','b'),axes):
    d=P[key]
    kw=dict(subsets=d['subsets'], set_labels=('','',''), set_colors=(GREEN,BLUE,ORANGE), alpha=0.42, ax=ax)
    if EQ: kw['layout_algorithm']=EQ
    v=venn3(**kw)
    ckw=dict(subsets=d['subsets'], linewidth=1.3, ax=ax)
    if EQ: ckw['layout_algorithm']=EQ
    circles=venn3_circles(**ckw)
    for c,col in zip(circles,(GREEN,BLUE,ORANGE)): c.set_edgecolor(col)
    for sid in ('100','010','001','110','101','011','111'):
        t=v.get_label_by_id(sid)
        if t:
            t.set_fontsize(17 if sid=='111' else 11.5)
            t.set_fontweight('bold' if sid=='111' else 'normal'); t.set_color('#111')
    ax.set_xlim(-1.08,1.08); ax.set_ylim(-1.18,1.28); ax.set_aspect('equal'); ax.axis('off')
    # set labels OUTSIDE circles, BELOW header
    ax.text(-0.54,0.66,'more isoforms',color=GREEN,fontsize=11.5,fontweight='bold',ha='center')
    ax.text(-0.54,0.575,f"({d['iso']})",color=GREEN,fontsize=9.5,ha='center')
    ax.text( 0.54,0.66,'more coding exons',color=BLUE,fontsize=12,fontweight='bold',ha='center')
    ax.text( 0.54,0.575,f"({d['exon']})",color=BLUE,fontsize=9.5,ha='center')
    ax.text( 0.0,-0.92,'longer CDS',color=ORANGE,fontsize=12,fontweight='bold',ha='center')
    ax.text( 0.0,-1.00,f"({d['cds']})",color=ORANGE,fontsize=9.5,ha='center')
    # rounded box + header + corner (axes coords)
    ax.add_patch(FancyBboxPatch((0.015,0.015),0.97,0.97, transform=ax.transAxes,
                 boxstyle='round,pad=0.004,rounding_size=0.025', fill=False, edgecolor='#9e9e9e', lw=1.1, zorder=0))
    ax.text(0.055,0.965,f"all {d['total']} revised {d['crop']} loci", transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top')
    ax.text(0.055,0.918,f"circles = {d['circ']} larger than the reference in ≥ 1 axis",
            transform=ax.transAxes, fontsize=9.5, color='#555', va='top')
    ax.text(0.05,0.11,d['other'], transform=ax.transAxes, fontsize=14, fontweight='bold', color='#333')
    ax.text(0.05,0.065,'revised in\nother ways', transform=ax.transAxes, fontsize=9, color='#555', va='top')
    ax.text(-0.02,1.10,key, transform=ax.transAxes, fontsize=15, fontweight='bold')
plt.subplots_adjust(wspace=0.06,left=0.02,right=0.98,top=0.95,bottom=0.02)
plt.savefig(OUT + '.png', dpi=150, bbox_inches='tight'); plt.savefig(OUT + '.svg', bbox_inches='tight')
print("done")
