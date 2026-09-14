#!/usr/bin/env python3
"""Figure 4a: GSE241782 bulk time-course module scores, publication style.
8 modules x 3 tissues; lines = tissue, x = time (glaucoma arm), mean +/- SEM.
Writes gse241782_module_tests.csv (delta, Welch 95% CI, nominal p, FDR)."""
import sys, os, gzip, re
sys.path.insert(0, '/mnt/agents/work/pipe/figs')
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import pubstyle as ps
ps.apply()

DATA = '/mnt/agents/work/pipe/data'
OUT = '/mnt/agents/work/pipe/out'
MODULES = {
 'Homeostatic_MG': ['P2ry12','Tmem119','Cx3cr1','Sall1','Siglech','Hexb','Fcrls','Olfml3','Gpr34','P2ry13'],
 'DAM_MG': ['Apoe','Spp1','Lgals3','Cst7','Trem2','Tyrobp','Gpnmb','Itgax','Cd9','Ctsb','Ctsd','Csf1'],
 'MHCII': ['Cd74','H2-Aa','H2-Ab1','H2-Eb1','Ciita'],
 'IFN_module': ['Ifit3','Isg15','Usp18','Ifit1','Oasl2','Irf7','Stat1'],
 'Pan_reactive': ['Gfap','Vim','Serpina3n','C3','C4b','Lcn2','Steap4','Osmr'],
 'A1_neurotoxic': ['H2-D1','H2-K1','Serping1','Ggta1','Iigp1','Gbp2','Fbln5','Ugt1a1','Fkbp5','Psmb8','Amigo2','Srgn'],
 'A2_neuroprotective': ['S100a10','Sphk1','Cd14','Ptx3','Tgm1','Emp1','Slc10a6','Tm4sf1','B3gnt5','Cd109','Ptgs2'],
 'Neutrophil': ['S100a8','S100a9','Ly6g','Mmp9','Cxcr2','Mpo'],
}
TISSUES = ['Retina','UON','MON']
TPS = ['0D','3D','2W','6W']

cnt = pd.read_csv(f'{DATA}/gse241782_counts.csv.gz')
gid = cnt.columns[0]
cnt[gid] = cnt[gid].astype(str).str.split('|').str[-1]
cnt = cnt.groupby(gid).sum()

# metadata: GSM -> title; column -> GSM via filelist
gsms=titles=None
with gzip.open(f'{DATA}/gse241782_matrix.txt.gz','rt') as f:
    for line in f:
        if line.startswith('!Sample_geo_accession'):
            gsms=[g.strip('"') for g in line.rstrip('\n').split('\t')[1:]]
        elif line.startswith('!Sample_title'):
            titles=[t.strip('"') for t in line.rstrip('\n').split('\t')[1:]]
meta = pd.DataFrame({'gsm':gsms,'title':titles})
parts = meta.title.str.split(',', expand=True)
meta['tissue']=parts[0].str.strip()
meta['tp']=parts[1].str.strip().replace({'Naïve':'0D'})
meta.loc[meta.title.str.contains('Naïve'),'tp']='0D'
meta['cond']=np.where(meta.title.str.contains('Glaucoma'),'OHT',
             np.where(meta.title.str.contains('Crush'),'Crush','Naive'))
# gsm -> column name from filelist
fl = pd.read_csv(f'{DATA}/gse241782_filelist.txt', sep='\t', comment='#', header=None,
                 names=['kind','name','date','size','type'])
fl = fl[fl.kind=='File']
fl['gsm']=fl['name'].str.extract(r'(GSM\d+)')
fl['col']=fl['name'].str.extract(r'GSM\d+_([A-Z]+-\d+)_genes')[0]
g2c = dict(zip(fl.gsm, fl.col))
meta['col']=meta.gsm.map(g2c)
assert meta.col.notna().all(), meta[meta.col.isna()]
print(meta.groupby(['tissue','tp','cond']).size())

lib = cnt.sum(axis=0)
lc = np.log2(cnt.div(lib,axis=1)*1e6 + 1)
S = {}
for mod, genes in MODULES.items():
    gs=[g for g in genes if g in lc.index]
    print(mod, len(gs),'/',len(genes))
    S[mod]=lc.loc[gs].mean(axis=0)
S = pd.DataFrame(S)

rows=[]
for _,r in meta.iterrows():
    for mod in MODULES:
        rows.append(dict(tissue=r.tissue, tp=r.tp, cond=r.cond, module=mod, score=S.loc[r.col, mod]))
df = pd.DataFrame(rows)

# glaucoma-arm tests vs tissue-matched naive (naive is shared across arms)
tres=[]
for tissue in TISSUES:
    base = df[(df.tissue==tissue)&(df.tp=='0D')].set_index('module')['score']
    for tp in ['3D','2W','6W']:
        for mod in MODULES:
            a = df[(df.tissue==tissue)&(df.tp==tp)&(df.cond=='OHT')&(df.module==mod)]['score']
            b = df[(df.tissue==tissue)&(df.tp=='0D')&(df.module==mod)]['score']
            t,p = stats.ttest_ind(a,b,equal_var=False)
            d = a.mean()-b.mean()
            na,nb=len(a),len(b)
            se = np.sqrt(a.var(ddof=1)/na + b.var(ddof=1)/nb)
            dfree = se**4/((a.var(ddof=1)/na)**2/(na-1)+(b.var(ddof=1)/nb)**2/(nb-1)) if se>0 else np.nan
            tcrit = stats.t.ppf(0.975, dfree) if dfree==dfree else np.nan
            tres.append(dict(tissue=tissue,tp=tp,module=mod,delta=d,lo=d-tcrit*se,hi=d+tcrit*se,p=p,n_case=na,n_ctrl=nb))
tres=pd.DataFrame(tres)
for tissue in TISSUES:
    m=tres.tissue==tissue
    tres.loc[m,'fdr']=stats.false_discovery_control(tres.loc[m,'p'])
tres.to_csv(f'{OUT}/gse241782_module_tests.csv',index=False)
print(tres[tres.p<0.1].sort_values('p').to_string())

# ---- figure: 2 rows x 4 cols, panel = module, lines = tissue
TC = {'Retina': ps.OI['vermillion'], 'UON': ps.OI['green'], 'MON': ps.OI['blue']}
fig, axes = plt.subplots(2,4, figsize=(7.2,3.4))
for ax,(mod,_) in zip(axes.ravel(), MODULES.items()):
    for tissue in TISSUES:
        sub = df[(df.tissue==tissue)&(df.module==mod)&((df.cond=='OHT')|(df.tp=='0D'))]
        means=[sub[sub.tp==tp]['score'].mean() for tp in TPS]
        sems=[sub[sub.tp==tp]['score'].sem() for tp in TPS]
        ax.errorbar(range(4), means, yerr=sems, color=TC[tissue], marker='o', ms=2.5,
                    lw=0.9, elinewidth=0.6, capsize=1.5, label=tissue)
        # stars for OHT vs naive
        for i,tp in enumerate(TPS[1:],1):
            r = tres[(tres.tissue==tissue)&(tres.tp==tp)&(tres.module==mod)].iloc[0]
            if r.p<0.05:
                ax.text(i, means[i]+ (sems[i] or 0)+ 0.02*(ax.get_ylim()[1]-ax.get_ylim()[0]),
                        ps.stars(r.p), ha='center', va='bottom', fontsize=5.5, color=TC[tissue])
    ax.set_xticks(range(4)); ax.set_xticklabels(['naive','3D','2W','6W'])
    ax.set_title(mod.replace('_',' '), pad=3)
    ps.despine(ax)
axes[0,0].set_ylabel('Module score (log2 CPM+1)')
axes[1,0].set_ylabel('Module score (log2 CPM+1)')
for ax in axes[1]: ax.set_xlabel('Time after OHT induction')
handles=[plt.Line2D([],[],color=TC[t],marker='o',ms=3,lw=0.9,label=t) for t in TISSUES]
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0,0.5))
fig.tight_layout(w_pad=1.2)
fig.savefig(f'{OUT}/fig4a_new.png')
print('saved fig4a_new.png')
