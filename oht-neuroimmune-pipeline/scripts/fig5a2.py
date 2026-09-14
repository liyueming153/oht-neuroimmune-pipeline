import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np, pandas as pd, sys
sys.path.insert(0,'/mnt/agents/work/pipe/figs'); import pubstyle as ps
R='/mnt/agents/work/oht-neuroimmune-pipeline/results/'
order=['Homeostatic MG','DAM-like MG','Iron+ MG (Fth1+)','Proliferating MG','MHCII+ Mac/APC',
       'Monocyte','Neutrophil','B cell','T cell','gdT cell','NK cell']
sh=pd.read_csv(R+'liana_Retina_Sham_1000perms.csv.gz')
mb=pd.read_csv(R+'liana_Retina_MB_1000perms.csv.gz')
def mat(d):
    sp=d[(d.ligand_complex=='Spp1')&(d.receptor_complex=='Cd44')]
    sp=sp[(sp.specificity_rank<=0.05)&(sp.cellphone_pvals<=0.05)]
    m=pd.DataFrame(np.nan,index=order,columns=order)
    for r in sp.itertuples():
        if r.source in order and r.target in order:
            m.loc[r.source,r.target]=r.lr_means
    return m
m_sh=mat(sh); m_mb=mat(mb)
vmax=float(np.nanmax(m_mb.values))
cmap=mcolors.LinearSegmentedColormap.from_list('w2v',['#FFFFFF','#F4C7B8',ps.OI['vermillion']])
cmap.set_bad('#F7F7F7')
fig,axes=plt.subplots(1,2,figsize=(176/25.4,88/25.4))
for ax,(m,ttl,note) in zip(axes,[(m_sh,'Retina \u2013 Sham','0 edges pass the specificity criterion\n(0 of 12,973 candidates)'),
                                  (m_mb,'Retina \u2013 OHT','18 of 33 edges pass the specificity criterion\n(24,842 candidates)')]):
    im=ax.imshow(m.values.astype(float),cmap=cmap,vmin=0,vmax=vmax,aspect='equal')
    ax.set_xticks(range(len(order))); ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order,rotation=90,fontsize=5.4)
    ax.set_yticklabels(order,fontsize=5.4)
    ax.set_xticks(np.arange(-.5,len(order),1),minor=True); ax.set_yticks(np.arange(-.5,len(order),1),minor=True)
    ax.grid(which='minor',color='white',lw=0.6); ax.tick_params(which='both',length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(ttl,fontsize=8)
    if np.isnan(m.values.astype(float)).all():
        ax.text(5,5,note,ha='center',va='center',fontsize=6.5,color='#4A4A4A')
    else:
        ax.set_xlabel('Receiver',fontsize=7)
        ax.text(5,-2.6,note,ha='center',va='top',fontsize=6)
axes[0].set_ylabel('Sender',fontsize=7)
cb=fig.colorbar(im,ax=axes,fraction=0.025,pad=0.02); cb.set_label('LR-means',fontsize=7); cb.ax.tick_params(labelsize=6)
cb.outline.set_visible(False)
fig.tight_layout(pad=0.4)
fig.savefig('fig5a_new.png',dpi=600,bbox_inches='tight'); plt.close(fig); print('fig5a v3 ok', m_mb.notna().sum().sum())
