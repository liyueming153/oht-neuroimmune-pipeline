#!/usr/bin/env python3
"""Figures 4b (Müller module bars, SHAM vs OHT) & 4c (volcano) from GSE306785.
FC = log2 ratio of means; p = Welch on log2(x+1). Publication style."""
import sys, os
sys.path.insert(0, '/mnt/agents/work/pipe/figs')
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import pubstyle as ps
ps.apply()

DATA='/mnt/agents/work/pipe/data'; OUT='/mnt/agents/work/pipe/out'
d = pd.read_csv(f'{DATA}/gse306785_macs.csv.gz', index_col=0)
print(d.shape, list(d.columns))
sham=[c for c in d.columns if 'SHAM' in c.upper()]
oht=[c for c in d.columns if 'OHT' in c.upper()]
print(sham, oht)
lx = np.log2(d+1)
log2fc = np.log2((d[oht].mean(axis=1)+1)/(d[sham].mean(axis=1)+1))
pvals = pd.Series([stats.ttest_ind(lx.loc[g,oht], lx.loc[g,sham], equal_var=False).pvalue
                   for g in d.index], index=d.index)
fdr = pd.Series(stats.false_discovery_control(pvals.fillna(1)), index=d.index)
res = pd.DataFrame({'log2fc':log2fc,'p':pvals,'fdr':fdr})
res.to_csv(f'{OUT}/gse306785_de.csv')
print(res.loc[['Cx3cl1','Spp1','Gfap']] if 'Cx3cl1' in d.index else 'gene check')
print('FDR<0.05:', (fdr<0.05).sum())

MODULES = {
 'Homeostatic_MG': ['P2ry12','Tmem119','Cx3cr1','Sall1','Siglech','Hexb','Fcrls','Olfml3','Gpr34','P2ry13'],
 'DAM_MG': ['Apoe','Spp1','Lgals3','Cst7','Trem2','Tyrobp','Gpnmb','Itgax','Cd9','Ctsb','Ctsd','Csf1'],
 'IFN_module': ['Ifit3','Isg15','Usp18','Ifit1','Oasl2','Irf7','Stat1'],
 'Pan_reactive': ['Gfap','Vim','Serpina3n','C3','C4b','Lcn2','Steap4','Osmr'],
 'Chemokine': ['Ccl5','Ccl3','Ccl4','Cxcl10','Ccl2','Ccl12'],
}
scores={}
for mod,genes in MODULES.items():
    gs=[g for g in genes if g in lx.index]
    scores[mod]=lx.loc[gs].mean(axis=0)
S=pd.DataFrame(scores)

# ---- Fig 4b: grouped bars SHAM vs OHT per module + per-sample dots
fig,ax=plt.subplots(figsize=(4.4,2.2))
x=np.arange(len(MODULES)); w=0.36
for i,(grp,cols,lab) in enumerate([( 'SHAM',sham,'SHAM'),('OHT',oht,'OHT')]):
    vals=[S.loc[cols,m].mean() for m in MODULES]
    errs=[S.loc[cols,m].sem() for m in MODULES]
    b=ax.bar(x+(i-0.5)*w, vals, w, yerr=errs, capsize=2,
             color=ps.COND['Sham' if grp=='SHAM' else 'OHT'],
             error_kw=dict(elinewidth=0.6), label=grp)
    for j,m in enumerate(MODULES):
        ax.scatter(np.full(len(cols), x[j]+(i-0.5)*w), S.loc[cols,m],
                   s=4, color='white', edgecolor='black', linewidth=0.4, zorder=3)
ax.set_xticks(x); ax.set_xticklabels([m.replace('_',' ') for m in MODULES], rotation=20, ha='right')
ax.set_ylabel('Module score (log2 norm. +1)')
ax.legend()
ps.despine(ax)
fig.tight_layout()
fig.savefig(f'{OUT}/fig4b_new.png'); print('saved fig4b')

# ---- Fig 4c: volcano
fig,ax=plt.subplots(figsize=(3.2,2.6))
neglog=-np.log10(res.p.clip(lower=1e-20))
ax.scatter(res.log2fc, neglog, s=2, color='#BBBBBB', linewidths=0, rasterized=True)
key=['Cx3cl1','Spp1','Gfap','Cx3cr1','Trem2','Apoe','Lgals3','Ccl5']
sig=(res.p<0.05)&(res.log2fc.abs()>1)
ax.scatter(res.loc[sig,'log2fc'], neglog[sig], s=3, color=ps.OI['vermillion'], linewidths=0, rasterized=True)
for g in key:
    if g in res.index:
        ax.annotate(g, (res.loc[g,'log2fc'], neglog[g]), fontsize=6, style='italic',
                    color=ps.OI['vermillion'] if g in ('Cx3cl1','Spp1') else 'black',
                    xytext=(3,3), textcoords='offset points')
ax.axhline(-np.log10(0.05), ls='--', lw=0.6, color='grey')
ax.axvline(0, lw=0.5, color='black')
ax.set_xlabel('log2 fold-change (OHT / SHAM)')
ax.set_ylabel('-log10 P (Welch)')
ps.despine(ax)
fig.tight_layout()
fig.savefig(f'{OUT}/fig4c_new.png'); print('saved fig4c')
