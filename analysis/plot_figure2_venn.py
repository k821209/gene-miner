#!/usr/bin/env python3
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
P = {
 'a': dict(crop='rice', total='20,717', circ='18,401', other='2,316',
           subsets=(203,382,639,742,240,1590,14605), iso='15,687', exon='17,216', cds='17,177'),
 'b': dict(crop='soybean', total='21,639', circ='9,213', other='12,426',
           subsets=(165,604,558,1404,98,2279,4105), iso='4,976', exon='7,546', cds='7,886'),
}
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
    ax.set_xlim(-1.0,1.0); ax.set_ylim(-1.15,1.25); ax.set_aspect('equal'); ax.axis('off')
    # set labels OUTSIDE circles, BELOW header
    ax.text(-0.66,0.62,'more isoforms',color=GREEN,fontsize=12,fontweight='bold',ha='center')
    ax.text(-0.66,0.53,f"({d['iso']})",color=GREEN,fontsize=9.5,ha='center')
    ax.text( 0.66,0.62,'more coding exons',color=BLUE,fontsize=12,fontweight='bold',ha='center')
    ax.text( 0.66,0.53,f"({d['exon']})",color=BLUE,fontsize=9.5,ha='center')
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
plt.savefig('/tmp/fig2_new.png',dpi=150,bbox_inches='tight'); plt.savefig('/tmp/fig2_new.svg',bbox_inches='tight')
print("done")
