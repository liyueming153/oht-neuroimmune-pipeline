import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, sys
from scipy.stats import wilcoxon, mannwhitneyu, ttest_ind
sys.path.insert(0,'/mnt/agents/work/pipe/figs'); import pubstyle as ps

feat=pd.read_csv('/mnt/agents/work/gse293358/GSM8880855_AMf20_features.tsv.gz',sep='\t',header=None,names=['id','name','type'])
emap=dict(zip(feat.id, feat.name))
d=pd.read_csv('/mnt/agents/work/pipe/data/gse304931_counts.csv.gz')
d['ens']=d['gene'].str.split('.').str[0]
d['sym']=d['ens'].map(emap)
d=d.dropna(subset=['sym'])
X=d.groupby('sym')[[c for c in d.columns if '_' in c]].mean()
X=np.log2(X+1)
WT=[f'WT_{i}' for i in range(1,5)]; IOP=[f'IOP_{i}' for i in range(1,5)]
KO=[f'Spp1KO_{i}' for i in range(1,5)]; KOI=[f'Spp1KO_IOP_{i}' for i in range(1,5)]
wt_fc=X[IOP].mean(1)-X[WT].mean(1)
ko_fc=X[KOI].mean(1)-X[KO].mean(1)
res=pd.DataFrame({'wt_fc':wt_fc,'ko_fc':ko_fc})

# key gene stats for text
for g in ['Spp1','Cx3cr1']:
    a=X.loc[g,IOP].values; b=X.loc[g,WT].values
    print(g,'WT+IOP log2FC=%.2f p=%.4f'%(a.mean()-b.mean(), ttest_ind(a,b,equal_var=False).pvalue))
a=X.loc['Cx3cr1',KOI].values; b=X.loc['Cx3cr1',IOP].values
print('Cx3cr1 KO+IOP vs WT+IOP log2FC=%.2f p=%.2e'%(a.mean()-b.mean(), ttest_ind(a,b,equal_var=False).pvalue))

mods={'Homeostatic':['P2ry12','Tmem119','Sall1','Cx3cr1','Hexb','Siglech'],
 'DAM':['Trem2','Apoe','Spp1','Gpnmb','Lgals3','Lpl','Cd9','Itgax','Clec7a','Cst7','Axl','Tyrobp'],
 'IFN':['Ifit1','Ifit2','Ifit3','Isg15','Irf7','Mx1','Mx2','Oas2','Oas3','Stat1','Rsad2'],
 'MHCII':['Cd74','H2-Aa','H2-Ab1','H2-Eb1','H2-Eb2','H2-DMb2','H2-DMa'],
 'Chemokine':['Ccl2','Ccl3','Ccl4','Ccl5','Ccl7','Ccl12','Cxcl10','Cxcl16']}
rows=[]
for m,gs in mods.items():
    gs=[g for g in gs if g in res.index]
    w=res.loc[gs,'wt_fc']; k=res.loc[gs,'ko_fc']
    p2=wilcoxon(k,w).pvalue
    p1=mannwhitneyu(k,w,alternative='greater').pvalue
    rows.append({'module':m,'n':len(gs),'wt_mean':w.mean(),'ko_mean':k.mean(),'wilcoxon_2sided':p2,'dependence_1sided':p1})
    print(m,'n=%d WT=%.2f KO=%.2f p2=%.4f p1=%.4f'%(len(gs),w.mean(),k.mean(),p2,p1))
df=pd.DataFrame(rows)
res.reset_index().rename(columns={'index':'sym'}).to_csv('gse304931_dependence_genes.csv',index=False)
df.to_csv('gse304931_dependence.csv',index=False)

# fig5c
fig,ax=plt.subplots(figsize=(96/25.4,76/25.4))
x=np.arange(len(df)); wd=0.36
ax.bar(x-wd/2,df.wt_mean,wd,color=ps.OI['blue'],label='WT+IOP vs WT')
ax.bar(x+wd/2,df.ko_mean,wd,color=ps.OI['vermillion'],label='$Spp1$-KO+IOP vs KO')
# p-value annotations removed: gene-level tests are descriptive (see Methods); statistics deposited in gse304931_dependence.csv
ax.axhline(0,color='k',lw=0.7)
ax.set_xticks(x); ax.set_xticklabels(df.module,fontsize=7)
ax.set_ylabel('Mean module log\u2082FC')
ax.set_ylim(top=df[['wt_mean','ko_mean']].max().max()*1.25)
ax.legend(frameon=False,fontsize=6,loc='upper left')
ps.despine(ax); fig.tight_layout(pad=0.3); fig.savefig('fig5c_new.png',dpi=600); plt.close(fig)

# fig5d: 8 genes x 4 groups
genes=['Spp1','Lgals3','Apoe','Gpnmb','Ifit3','Isg15','Cd74','Ccl5']
groups=[('WT',WT,'#8C8C8C'),('WT+IOP',IOP,ps.OI['blue']),('$Spp1$-KO',KO,'#F0B3A1'),('KO+IOP',KOI,ps.OI['vermillion'])]
fig,axes=plt.subplots(2,4,figsize=(176/25.4,88/25.4))
for ax,g in zip(axes.flat,genes):
    for j,(lab,cols,c) in enumerate(groups):
        v=X.loc[g,cols].values if g in X.index else np.zeros(len(cols))
        ax.bar(j,v.mean(),0.66,color=c,lw=0,yerr=v.std()/np.sqrt(len(v)),error_kw=dict(lw=0.6,capsize=1.5))
        ax.scatter(np.full(len(v),j)+np.linspace(-0.16,0.16,len(v)),v,s=4,color='k',lw=0,zorder=3)
    ax.set_title(g,fontsize=7.5,style='italic')
    ax.set_xticks(range(4)); ax.set_xticklabels(['WT','WT\n+IOP','KO','KO\n+IOP'],fontsize=5.6)
    ax.tick_params(axis='y',labelsize=5.6)
    ps.despine(ax)
axes[0,0].set_ylabel('log\u2082(norm+1)',fontsize=6.5); axes[1,0].set_ylabel('log\u2082(norm+1)',fontsize=6.5)
fig.tight_layout(pad=0.4,h_pad=0.8); fig.savefig('fig5d_new.png',dpi=600); plt.close(fig); print('fig5c/5d ok')
