#!/usr/bin/env python3
"""Step 05: Bulk time-course module replication (GSE241782; Methods 4.5).

Counts log-CPM transformed; OHT 3d/2W/6W vs naive contrasts (Welch t).
Module scores = means of log-normalized expression of module genes
(Supplementary Table S1). Sample columns mapped via the series SOFT file
and verified by Xist expression against recorded sex.
"""
import gzip, io
import os
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.environ.get('GSE241782_DIR', 'data/GSE241782')
OUT = os.environ.get('OUT_DIR', 'out')

MODULES = {
 'Homeostatic_MG': ['P2ry12','Tmem119','Cx3cr1','Sall1','Siglech','Hexb','Fcrls','Olfml3','Gpr34','P2ry13'],
 'DAM_MG': ['Apoe','Spp1','Lgals3','Cst7','Trem2','Tyrobp','Gpnmb','Itgax','Cd9','Ctsb','Ctsd','Csf1'],
 'MHCII': ['Cd74','H2-Aa','H2-Ab1','H2-Eb1','Ciita'],
 'IFN_module': ['Ifit3','Isg15','Usp18','Ifit1','Oasl2','Irf7','Stat1'],
 'Chemokine': ['Ccl5','Ccl3','Ccl4','Cxcl10','Ccl2','Ccl12'],
 'Iron+ MG': ['Fth1','Ftl1','Ltf','Hmox1','Slc40a1'],
 'T_NK_infil': ['Cd3e','Cd3g','Trbc1','Trbc2','Nkg7','Gzma','Klrb1c'],
 'Neutrophil': ['S100a8','S100a9','Ly6g','Mmp9','Cxcr2','Mpo'],
}

# counts: one CSV/MTX per GEO deposit; adapt loader to the deposited format
counts = pd.read_csv(f'{DATA}/GSE241782_counts.csv.gz', index_col=0)
meta = pd.read_csv(f'{DATA}/GSE241782_sample_meta.csv')  # column, tissue, condition, timepoint

def logcpm(df):
    lib = df.sum(axis=0)
    return np.log2(df.div(lib, axis=1) * 1e6 + 1)

lc = logcpm(counts)
scores = {}
for mod, genes in MODULES.items():
    gs = [g for g in genes if g in lc.index]
    scores[mod] = lc.loc[gs].mean(axis=0)
S = pd.DataFrame(scores)
S.to_csv(f'{OUT}/gse241782_module_scores.csv')

rows = []
for tissue in ['Retina', 'UON', 'MON']:
    for tp in ['3D', '2W', '6W']:
        case = meta[(meta.tissue == tissue) & (meta.timepoint == tp) & (meta.condition == 'OHT')]['column']
        ctrl = meta[(meta.tissue == tissue) & (meta.timepoint == 'naive')]['column']
        for mod in MODULES:
            a, b = S.loc[case, mod], S.loc[ctrl, mod]
            t, p = stats.ttest_ind(a, b, equal_var=False)
            rows.append(dict(tissue=tissue, timepoint=tp, module=mod,
                             delta=a.mean() - b.mean(), p=p))
res = pd.DataFrame(rows)
res['fdr'] = np.nan
for tissue in res.tissue.unique():
    m = res.tissue == tissue
    res.loc[m, 'fdr'] = stats.false_discovery_control(res.loc[m, 'p'])
res.to_csv(f'{OUT}/gse241782_module_tests.csv', index=False)
print(res.head())
