#!/usr/bin/env python3
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
# ---- data ----
# genome: reads, GM_rep, GM_alliso, BRAKER3
G = {
 'soybean':      dict(reads=0.89e8, rep=91.9, allo=92.1, brk=77.6, italic=False),
 'C. elegans':   dict(reads=1.12e8, rep=95.1, allo=98.3, brk=97.7, italic=True),
 'rice':         dict(reads=1.18e8, rep=93.6, allo=96.5, brk=90.2, italic=False),
 'D. melanogaster': dict(reads=1.39e8, rep=96.5, allo=98.9, brk=95.6, italic=True),
}
refs = [('rice RefSeq',8.87e9,99.2), ('soybean NCBI',1.069e10,99.3)]
GREEN='#2e7d32'; RED='#c0392b'; BLUE='#1f6fb2'
fig,(axa,axb)=plt.subplots(1,2,figsize=(13,5.4),gridspec_kw={'width_ratios':[1.05,1]})
# ===== Panel a =====
for name,d in G.items():
    x=d['reads']
    # dashed connector rep<->all-isoform
    axa.plot([x,x],[d['rep'],d['allo']],ls=':',color=GREEN,lw=1.2,zorder=1)
    axa.plot(x,d['rep'],'o',color=GREEN,ms=9,zorder=4)
    axa.plot(x,d['allo'],'o',mfc='none',mec=GREEN,mew=1.8,ms=9,zorder=4)
    axa.plot(x,d['brk'],'v',mfc='none',mec=RED,mew=1.8,ms=10,zorder=4)
for nm,x,y in refs:
    axa.plot(x,y,'D',color=BLUE,ms=12,zorder=4)
axa.axhline(99.0,ls=':',color=BLUE,lw=1,zorder=0)
axa.text(2.2e8,100.4,'≈ 99% curated-reference ceiling',color=BLUE,fontsize=8.5,ha='left')
# genome labels: leader lines to the GM-representative point, spread to the LEFT
lab_x=4.3e7
lab_y={'D. melanogaster':98.2,'C. elegans':95.7,'rice':93.4,'soybean':91.4}
for name,d in G.items():
    ly=lab_y[name]
    axa.annotate(name, xy=(d['reads'],d['rep']), xytext=(lab_x,ly),
                 fontsize=10, fontweight='bold', color=GREEN, va='center', ha='left',
                 style=('italic' if d['italic'] else 'normal'),
                 arrowprops=dict(arrowstyle='-',color='#9e9e9e',lw=0.8,shrinkA=2,shrinkB=6))
# reference callout
axa.annotate('rice & soybean\ncurated references', xy=(refs[0][1],99.2), xytext=(1.3e9,95.5),
             fontsize=9, color=BLUE, ha='center',
             arrowprops=dict(arrowstyle='-',color=BLUE,lw=0.9))
# evidence arrow
axa.annotate('', xy=(1.5e8,80.2), xytext=(7.5e9,80.2),
             arrowprops=dict(arrowstyle='<->',color='#333',lw=1.3))
axa.text(1.2e9,81.0,'~100× less evidence',fontsize=9.5,ha='center',color='#333')
axa.set_xscale('log'); axa.set_xlim(3.5e7,1.6e10); axa.set_ylim(75,101)
axa.set_xlabel('RNA-seq reads used for annotation (log scale)',fontsize=10)
axa.set_ylabel('BUSCO completeness (% complete)',fontsize=10)
axa.legend(handles=[
    Line2D([0],[0],marker='o',color='w',mfc=GREEN,ms=9,label='Gene-Miner, representative'),
    Line2D([0],[0],marker='o',color='w',mfc='none',mec=GREEN,mew=1.8,ms=9,label='Gene-Miner, all-isoform'),
    Line2D([0],[0],marker='v',color='w',mfc='none',mec=RED,mew=1.8,ms=10,label='BRAKER3 (same inputs)'),
    Line2D([0],[0],marker='D',color='w',mfc=BLUE,ms=11,label='Curated reference')],
    fontsize=8.5, loc='lower right', frameon=True, framealpha=0.95)
axa.set_title('a',loc='left',fontweight='bold',fontsize=13)
for s in ('top','right'): axa.spines[s].set_visible(False)
# ===== Panel b =====
genomes=['rice','soybean','D. melanogaster','C. elegans']
gm=[1.89,1.24,1.69,1.74]; brk=[1.20,1.14,1.28,1.32]; ref=[1.19,1.58,2.21,1.40]
x=np.arange(len(genomes)); w=0.26
b1=axb.bar(x-w,gm,w,color=GREEN,label='Gene-Miner')
b2=axb.bar(x,brk,w,color=RED,label='BRAKER3 (same reads)')
b3=axb.bar(x+w,ref,w,color=BLUE,label='Curated reference')
for bars,vals in ((b1,gm),(b2,brk),(b3,ref)):
    for bar,v in zip(bars,vals):
        axb.text(bar.get_x()+bar.get_width()/2, v+0.03, f'{v:.2f}',ha='center',fontsize=8)
axb.axhline(1.0,ls=':',color='#888',lw=1)
axb.text(3.48,1.03,'1 isoform / gene',fontsize=8.5,color='#666',ha='right',va='bottom')
axb.set_xticks(x); axb.set_xticklabels([g if ' ' not in g else '$\\it{%s}$'%g.replace(' ','\\ ') for g in genomes],fontsize=9)
axb.set_ylabel('Transcripts per gene',fontsize=10); axb.set_ylim(0,2.5)
axb.legend(fontsize=8.5,loc='upper left',frameon=True)
axb.set_title('b',loc='left',fontweight='bold',fontsize=13)
for s in ('top','right'): axb.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig('/tmp/fig4_new.png',dpi=150); plt.savefig('/tmp/fig4_new.svg')
print("done")
