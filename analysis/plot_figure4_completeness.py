#!/usr/bin/env python3
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import numpy as np
# Every plotted number is read from the analysis outputs, not typed in:
#   BUSCO from each runs short_summary, transcripts per gene from gm_compare.pys
#   report (catalogue and reference) and from the same report run on BRAKER3.
import glob, os, re
R = os.path.expanduser("~/gene-miner-runs"); REV = R + "/revision"
def busco(d):
    f = glob.glob(os.path.join(d, "short_summary*.txt"))[0]
    return float(re.search(r"C:([0-9.]+)%", open(f).read()).group(1))
def gmc(path, key):
    for line in open(path):
        k, *v = line.rstrip("\n").split("\t")
        if k.strip() == key: return float(v[0])
SRC = {  # name: (catalogue dir, BRAKER3 BUSCO dir, read count, italic)
 "soybean":         (REV + "/regen/soybean", R + "/soybean/braker_soy/busco_braker_soy_primary", 1.785e8, False),
 "C. elegans":      (REV + "/regen/cele",    R + "/cele/braker_cele/busco_braker_primary",       2.246e8, True),
 "rice":            (REV + "/regen/rice_s",  R + "/rice/etp_run/busco_braker_full",              2.363e8, False),
 "D. melanogaster": (REV + "/regen/dmel",    R + "/dmel/braker_dmel/busco_braker_dmel_primary",  2.784e8, True),
}
BRK_TSV = {"soybean": "soybean", "C. elegans": "cele", "rice": "rice", "D. melanogaster": "dmel"}
G = {}
for name, (cat, bbusco, reads, it) in SRC.items():
    G[name] = dict(reads=reads, rep=busco(cat + "/busco_primary"), allo=busco(cat + "/busco_alliso"),
                   brk=busco(bbusco), italic=it,
                   tpg=gmc(cat + "/gm_compare.tsv", "gene_miner_isoforms_per_gene"),
                   ref_tpg=gmc(cat + "/gm_compare.tsv", "reference_isoforms_per_gene"),
                   brk_tpg=gmc(REV + "/braker/" + BRK_TSV[name] + ".tsv", "gene_miner_isoforms_per_gene"))
    print(name, G[name])
refs = [('rice RefSeq',8.87e9,99.2), ('soybean NCBI',1.069e10,99.3)]
GREEN='#2e7d32'; RED='#c0392b'; BLUE='#1f6fb2'
fig=plt.figure(figsize=(13.4,5.8))
outer=gridspec.GridSpec(1,2,width_ratios=[1.18,1.0],wspace=0.22)
la=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=outer[0],width_ratios=[3.2,1.0],wspace=0.06)
axL=fig.add_subplot(la[0]); axR=fig.add_subplot(la[1],sharey=axL)
axb=fig.add_subplot(outer[1])
def draw_genomes(ax):
    for name,d in G.items():
        x=d['reads']
        ax.axvline(x,color='#cfcfcf',lw=0.9,zorder=0)
        ax.plot([x,x],[d['rep'],d['allo']],ls=':',color=GREEN,lw=1.2,zorder=2)
        ax.plot(x,d['rep'],'o',color=GREEN,ms=9,zorder=4)
        ax.plot(x,d['allo'],'o',mfc='none',mec=GREEN,mew=1.8,ms=9,zorder=4)
        ax.plot(x,d['brk'],'v',mfc='none',mec=RED,mew=1.8,ms=10,zorder=4)
draw_genomes(axL)
for nm,x,y in refs:
    axR.axvline(x,color='#dfe7ef',lw=0.9,zorder=0)
    axR.plot(x,y,'D',color=BLUE,ms=12,zorder=4)
for ax in (axL,axR):
    ax.axhline(99.0,ls=':',color=BLUE,lw=1,zorder=0)
    ax.set_ylim(75,101)
# left = the four genomes, linear zoom so they separate
axL.set_xlim(1.60e8,3.00e8)
axL.set_xticks([]);
transL=blended_transform_factory(axL.transData,axL.transAxes)
def readfmt(r):
    return f'{r/1e9:.2f}B' if r>=1e9 else f'{r/1e6:.0f}M'
def lbl(name,r,italic):
    nm = ('$\\it{%s}$'%name.replace(' ','\\ ')) if italic else name
    return f'{nm} ({readfmt(r)})'
for name,d in G.items():
    axL.text(d['reads'],-0.03,lbl(name,d['reads'],d['italic']),transform=transL,rotation=30,
             ha='right',va='top',fontsize=9.2,color=GREEN)
# right = curated references
axR.set_xlim(7.0e9,1.30e10)
axR.set_xticks([])
transR=blended_transform_factory(axR.transData,axR.transAxes)
for nm,x,y in refs:
    axR.text(x,-0.03,f'{nm} ({readfmt(x)})',transform=transR,rotation=30,ha='right',va='top',
             fontsize=8.8,color=BLUE)
# hide inner spines and draw axis-break marks
axL.spines['right'].set_visible(False); axR.spines['left'].set_visible(False)
axR.tick_params(labelleft=False,left=False)
for ax in (axL,axR):
    ax.spines['top'].set_visible(False)
d=.012
kw=dict(transform=axL.transAxes,color='k',clip_on=False,lw=1.1)
axL.plot([1-d,1+d],[-d,d],**kw); axL.plot([1-d,1+d],[1-d,1+d],**kw)
kw=dict(transform=axR.transAxes,color='k',clip_on=False,lw=1.1)
axR.plot([-d*3.2,d*3.2],[-d,d],**kw); axR.plot([-d*3.2,d*3.2],[1-d,1+d],**kw)
# ceiling label
axL.text(1.63e8,100.3,'≈ 99% curated-reference ceiling',color=BLUE,fontsize=8.5,ha='left')
axL.set_ylabel('BUSCO completeness (% complete)',fontsize=10)
axL.set_title('a',loc='left',fontweight='bold',fontsize=13)
axL.legend(handles=[
    Line2D([0],[0],marker='o',color='w',mfc=GREEN,ms=9,label='Gene-Miner, representative'),
    Line2D([0],[0],marker='o',color='w',mfc='none',mec=GREEN,mew=1.8,ms=9,label='Gene-Miner, all-isoform'),
    Line2D([0],[0],marker='v',color='w',mfc='none',mec=RED,mew=1.8,ms=10,label='BRAKER3 (same inputs)'),
    Line2D([0],[0],marker='D',color='w',mfc=BLUE,ms=11,label='Curated reference')],
    fontsize=8.3,loc='lower left',bbox_to_anchor=(0.26,0.015),frameon=True,framealpha=0.95)
fig.text(0.30,0.002,'RNA-seq reads used for annotation',fontsize=10,ha='center')
for s in ('top','right'): axb.spines[s].set_visible(False)
# ===== Panel b =====
genomes=['rice','soybean','D. melanogaster','C. elegans']
gm=[G[g]['tpg'] for g in genomes]; brk=[G[g]['brk_tpg'] for g in genomes]; ref=[G[g]['ref_tpg'] for g in genomes]
x=np.arange(len(genomes)); w=0.26
b1=axb.bar(x-w,gm,w,color=GREEN,label='Gene-Miner')
b2=axb.bar(x,brk,w,color=RED,label='BRAKER3 (same reads)')
b3=axb.bar(x+w,ref,w,color=BLUE,label='Curated reference')
for bars,vals in ((b1,gm),(b2,brk),(b3,ref)):
    for bar,v in zip(bars,vals):
        axb.text(bar.get_x()+bar.get_width()/2,v+0.03,f'{v:.2f}',ha='center',fontsize=8)
axb.axhline(1.0,ls=':',color='#888',lw=1)
axb.set_xticks(x); axb.set_xticklabels([g if ' ' not in g else '$\\it{%s}$'%g.replace(' ','\\ ') for g in genomes],fontsize=9)
axb.set_ylabel('Transcripts per gene',fontsize=10); axb.set_ylim(0,2.5)
h, l = axb.get_legend_handles_labels()
h.append(Line2D([0],[0],ls=':',color='#888',lw=1)); l.append('one transcript per gene')
axb.legend(h, l, fontsize=8.5,loc='upper left',frameon=True)
axb.set_title('b',loc='left',fontweight='bold',fontsize=13)
plt.subplots_adjust(left=0.06,right=0.985,top=0.94,bottom=0.22)
# clean '~100x more evidence' arrow spanning the break, in figure coords
posL=axL.get_position(); posR=axR.get_position()
ax0=posL.x0+0.62*(posL.x1-posL.x0)   # over the genome cluster
ax1=posR.x0+0.55*(posR.x1-posR.x0)   # over the curated refs
ay=posL.y0+0.40*(posL.y1-posL.y0)
fig.add_artist(plt.matplotlib.patches.FancyArrowPatch((ax0,ay),(ax1,ay),
    transform=fig.transFigure,arrowstyle='<->',mutation_scale=14,color='#333',lw=1.3))
fig.text((ax0+ax1)/2,ay+0.02,'~40–60× more sequencing evidence',ha='center',fontsize=9.3,color='#333')
import sys
out = sys.argv[1] if len(sys.argv) > 1 else '/tmp/fig4_new.png'
plt.savefig(out,dpi=160); print("done",out)
