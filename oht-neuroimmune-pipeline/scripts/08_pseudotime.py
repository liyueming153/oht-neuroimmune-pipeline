#!/usr/bin/env python3
"""Step 08: Microglial sub-clustering and diffusion pseudotime (Methods 4.3).

Six microglial classes re-clustered (Leiden 0.8); states scored against
eight modules (homeostatic, DAM, MHCII, interferon, chemokine, iron,
proliferating, ribosomal). Diffusion maps (12 components) and DPT rooted
at the maximally homeostatic cell. Trajectory genes by Spearman
correlation with DPT (8,000-cell subsample).
"""
import numpy as np
import os
import pandas as pd
import scanpy as sc

OUT = os.environ.get('OUT_DIR', 'out')
adata = sc.read_h5ad(f'{OUT}/annotated.h5ad')

MG = ['Homeostatic MG', 'DAM-like MG', 'Iron+ MG (Fth1+)',
      'Activated MG (Ccl5+)', 'Proliferating MG', 'BAM']
mg = adata[adata.obs.cell_class.isin(MG)].copy()
print('microglial cells:', mg.n_obs)

sc.pp.highly_variable_genes(mg, n_top_genes=2500, flavor='seurat', batch_key='sample')
sc.pp.scale(mg, zero_center=False, max_value=10)
sc.pp.pca(mg, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(mg, n_neighbors=15, n_pcs=30)
sc.tl.umap(mg)
sc.tl.leiden(mg, resolution=0.8, key_added='mg_leiden')

MODULES = {
 'homeostatic': ['P2ry12','Tmem119','Cx3cr1','Sall1','Siglech','Hexb','Fcrls','Olfml3','Gpr34','P2ry13'],
 'DAM': ['Apoe','Spp1','Lgals3','Cst7','Trem2','Tyrobp','Gpnmb','Itgax','Cd9','Ctsb','Ctsd','Csf1'],
 'MHCII': ['Cd74','H2-Aa','H2-Ab1','H2-Eb1','Ciita'],
 'interferon': ['Ifit3','Isg15','Usp18','Ifit1','Oasl2','Irf7','Stat1'],
 'chemokine': ['Ccl5','Ccl3','Ccl4','Cxcl10','Ccl2','Ccl12'],
 'iron': ['Fth1','Ftl1','Ltf','Hmox1','Slc40a1'],
 'proliferating': ['Mki67','Top2a','Birc5','Cdk1'],
 'ribosomal': ['Rpl3','Rplp0','Rps3','Rps18','Rps27'],
}
for name, genes in MODULES.items():
    gs = [g for g in genes if g in mg.var_names]
    sc.tl.score_genes(mg, gs, score_name=f'mod_{name}', use_raw=False)

# diffusion pseudotime rooted at the maximally homeostatic cell
sc.tl.diffmap(mg, n_comps=12)
mg.uns['iroot'] = int(np.argmax(mg.obs['mod_homeostatic'].values))
sc.tl.dpt(mg, n_dcs=12)

# trajectory genes: Spearman correlation with DPT on 8,000-cell subsample
rng = np.random.default_rng(1337)
take = rng.choice(mg.n_obs, size=min(8000, mg.n_obs), replace=False)
sub = mg[take]
dpt = sub.obs['dpt_pseudotime'].values
X = sub.X.toarray() if hasattr(sub.X, 'toarray') else sub.X
from scipy.stats import spearmanr
rhos = np.array([spearmanr(X[:, j], dpt).statistic for j in range(X.shape[1])])
traj = pd.DataFrame({'gene': sub.var_names, 'rho': rhos}).dropna()
traj['abs_rho'] = traj.rho.abs()
traj.sort_values('abs_rho', ascending=False).to_csv(f'{OUT}/dpt_trajectory_genes.csv', index=False)

mg.write(f'{OUT}/microglia_dpt.h5ad', compression='gzip')
print(traj.sort_values('abs_rho', ascending=False).head(20))
