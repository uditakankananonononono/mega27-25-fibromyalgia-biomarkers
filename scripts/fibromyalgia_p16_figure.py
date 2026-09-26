"""Render frozen P16 direction reversal with post-outcome sex sensitivity, no biomarker claim."""
import json
import numpy as np
import matplotlib.pyplot as plt
x=json.load(open('results/fibromyalgia_p16_result.json'))
a=[z for z in x['frozen'] if z['measured']];sx={z['gene']:z for z in x['post_outcome_sex_sensitivity']}
names=[z['gene'] for z in a];y=np.arange(len(names))[::-1]
fig,ax=plt.subplots(figsize=(7.6,4.4))
ax.axvline(0,color='#555',linewidth=.8)
for i,z in enumerate(a):
 row=y[i];ax.plot(z['normal_95_CI'],[row+.12]*2,color='#205b65',linewidth=1.8)
 ax.scatter(z['g'],row+.12,color='#205b65',s=35,zorder=4)
 ax.scatter(sx[z['gene']]['female_g'],row-.12,color='#bd7049',marker='D',s=26,zorder=4)
ax.set(yticks=y,yticklabels=names,xlabel='Hedges g (fibromyalgia minus control)',xlim=(-.5,1.13))
ax.text(-.46,max(y)+.58,'Frozen GSE67311 direction: DOWN',color='#555',fontsize=9)
ax.scatter([],[],color='#205b65',label='All 96 FM / 93 controls (95% normal interval)')
ax.scatter([],[],color='#bd7049',marker='D',label='Female-only 91 FM / 41 controls (post-outcome)')
ax.legend(frameon=False,loc='lower right',fontsize=8)
ax.set_ylim(-.7,max(y)+1.0)
ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig('paper/figures/fibromyalgia_p16_effects.pdf');plt.close(fig)
