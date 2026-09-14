#!/usr/bin/env python3
"""Step 02c: final hierarchical annotation (v3).

Base: per-cluster argmax of sc.tl.score_genes z-scores over curated
cell-class signatures (as in Methods 4.2).
Overrides via per-cluster raw mean log-normalized expression:
contaminant gate; lymphocyte/myeloid lineage gates.
Cell-level refinement splits mixed lymphocyte clusters by Cd3e/Trdc/Nkg7.
"""
import os, gc
import numpy as np
import pandas as pd
import scanpy as sc

OUT = os.environ.get('OUT_DIR', 'out')
_pf = f'{OUT}/clustered_ckpt.h5ad'; _ps = f'{OUT}/clustered_slim.h5ad'
adata = sc.read_h5ad(_ps if os.path.exists(_ps) else _pf)  # prefer slim
cl = adata.obs['leiden1']

ZSIG = {
 'Homeostatic MG':       ['P2ry12','Tmem119','Cx3cr1','Sall1','Siglech','Hexb','Fcrls','Olfml3','Gpr34','P2ry13'],
 'DAM-like MG':          ['Apoe','Spp1','Lgals3','Cst7','Trem2','Tyrobp','Gpnmb','Itgax','Cd9','Ctsb'],
 'Iron+ MG (Fth1+)':     ['Fth1','Ftl1','Ltf','Hmox1','Slc40a1'],
 'Activated MG (Ccl5+)': ['Ccl5','Ccl3','Ccl4','Cxcl10','Ccl2','Ccl12'],
 'Proliferating MG':     ['Mki67','Top2a','Birc5','Cdk1'],
 'BAM':                  ['Mrc1','Lyve1','Pf4','Folr2','Cd163'],
 'MHCII+ Mac/APC':       ['Cd74','H2-Aa','H2-Ab1','H2-Eb1','Ciita'],
 'Monocyte':             ['Ly6c2','Ccr2','Plac8','Chil3'],
 'pDC':                  ['Bst2','Cox6a2','Klk1','Tcf4'],
 'T cell':               ['Cd3e','Cd3g','Trbc1','Trbc2'],
 'gdT cell':             ['Trdc','Trgc1','Trgc2'],
 'NK cell':              ['Nkg7','Klrb1c','Ncr1','Gzma'],
 'B cell':               ['Cd79a','Cd79b','Ms4a1','Cd19'],
 'Neutrophil':           ['S100a8','S100a9','Ly6g','Cxcr2','Mmp9','Mpo'],
}
RAW = {  # for gates
 'B': ['Cd79a','Cd79b','Ms4a1','Cd19'], 'gdT': ['Trdc','Trgc1','Trgc2'],
 'T': ['Cd3e','Cd3g','Trbc1','Trbc2'], 'NK': ['Nkg7','Klrb1c','Ncr1','Gzma'],
 'Neut': ['S100a8','S100a9','Ly6g','Cxcr2','Mmp9'],
 'pDC': ['Bst2','Cox6a2','Klk1','Tcf4'], 'Mono': ['Ly6c2','Ccr2','Plac8','Chil3'],
 'HomeoMG': ['P2ry12','Tmem119','Cx3cr1','Sall1','Hexb'],
}
CONTAM = {
 'rod/cone': ['Rho','Pde6a','Pde6b','Nrl','Opn1sw'], 'RGC': ['Rbpms','Nefl','Sncg','Slc17a6'],
 'Muller': ['Rlbp1','Glul','Slc1a3','Clu'], 'astrocyte': ['Gfap','Aqp4','S100b','Aldh1l1'],
 'endothelial': ['Pecam1','Kdr','Cldn5'], 'RPE': ['Rpe65','Mitf','Tyr'],
}

def cl_mean_raw(genes):
    gs = [g for g in genes if g in adata.var_names]
    X = adata[:, gs].X
    df = pd.DataFrame(X.toarray() if hasattr(X, 'toarray') else X, columns=gs)
    df['cl'] = cl.values
    return df.groupby('cl', observed=True).mean().mean(axis=1)

print('z-scoring signatures...', flush=True)
zsc = {}
for name, genes in ZSIG.items():
    gs = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, gs, score_name=f'z_{name}', use_raw=False)
    zsc[name] = adata.obs[f'z_{name}']
Z = pd.DataFrame(zsc)
Zc = Z.groupby(cl.values).mean()

raw = pd.DataFrame({k: cl_mean_raw(v) for k, v in RAW.items()})
con = pd.DataFrame({k: cl_mean_raw(v) for k, v in CONTAM.items()})

assign = {}
for c in Zc.index:
    immune_raw = raw.loc[c].max()
    if immune_raw < 0.35 and con.loc[c].max() > immune_raw:
        assign[c] = 'Contaminant'; continue
    R = raw.loc[c]
    if R['B'] > 0.5:   assign[c] = 'B cell'; continue
    if R['gdT'] > 0.5: assign[c] = 'gdT-mixed'; continue
    if R['T'] > 0.5 and R['NK'] > 0.4: assign[c] = 'TNK-mixed'; continue
    if R['T'] > 0.5:   assign[c] = 'T-mixed'; continue
    if R['NK'] > 0.5 and R['T'] < 0.3: assign[c] = 'NK cell'; continue
    if R['Neut'] > 0.5: assign[c] = 'Neutrophil'; continue
    if R['pDC'] > 0.5 and R['HomeoMG'] < 0.8: assign[c] = 'pDC'; continue
    if R['Mono'] > 0.5 and R['HomeoMG'] < 1.0: assign[c] = 'Monocyte'; continue
    assign[c] = Zc.loc[c].idxmax()

print('cluster assignments:'); print(pd.Series(assign).value_counts(), flush=True)
adata.obs['cell_class'] = cl.map(assign).astype(str)

# cell-level refinement of lymphocyte clusters
mixed = [c for c, v in assign.items() if v in ('gdT-mixed', 'TNK-mixed', 'T-mixed')]
if mixed:
    m = adata.obs['leiden1'].isin(mixed).values
    sub = adata[m]
    def expr(g):
        return sub[:, g].X.toarray().ravel() if g in sub.var_names else np.zeros(sub.n_obs)
    cd3, trdc, nkg7 = expr('Cd3e'), expr('Trdc'), expr('Nkg7')
    new = np.where(trdc > 0, 'gdT cell', np.where(cd3 > 0, 'T cell', 'NK cell'))
    cur = adata.obs['cell_class'].astype(str)
    cur.loc[adata.obs.index[m]] = new
    adata.obs['cell_class'] = cur
    print('refined', int(m.sum()), 'lymphoid cells:', pd.Series(new).value_counts().to_dict(), flush=True)

adata.obs['cell_class'] = adata.obs['cell_class'].astype('category')
adata.obs['is_contaminant'] = (adata.obs['cell_class'] == 'Contaminant')
dropcols = [c for c in adata.obs.columns if c.startswith(('sc_', 'z_'))]
adata.obs.drop(columns=dropcols, inplace=True, errors='ignore')
Z.to_csv(f'{OUT}/cell_class_zscores.csv.gz')
adata.write(f'{OUT}/annotated.h5ad', compression='lzf')
print('FINAL (excl. contaminants):', int((~adata.obs.is_contaminant).sum()), flush=True)
print(pd.crosstab(adata.obs.cell_class, [adata.obs.tissue, adata.obs.condition]), flush=True)
