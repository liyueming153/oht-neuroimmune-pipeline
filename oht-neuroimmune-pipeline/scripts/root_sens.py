#!/usr/bin/env python3
"""Root-sensitivity analysis for DPT (promised in Methods 4.3).
Re-derive DPT from the 10 most-homeostatic cells individually and from
5 random homeostatic roots; report Spearman rho vs the primary DPT."""
import sys, os
import numpy as np, pandas as pd, scanpy as sc
OUT = os.environ.get('OUT_DIR', '/tmp/ohtrun/out')
mg = sc.read_h5ad(f'{OUT}/microglia_dpt.h5ad')
base = mg.obs['dpt_pseudotime'].values.copy()
homeo_score = mg.obs['mod_homeostatic'].values
order = np.argsort(-homeo_score)
rng = np.random.default_rng(42)
top10 = order[:10]
# random roots from homeostatic-labeled cells
hmask = (mg.obs['mg_leiden'].map(lambda c: True).values)  # placeholder
# homeostatic cells: mod_homeostatic above median among cells whose module max is homeostatic
mods = ['homeostatic','DAM','MHCII','interferon','chemokine','iron','proliferating','ribosomal']
Z = mg.obs[[f'mod_{m}' for m in mods]].values
is_h = Z.argmax(axis=1)==0
hcells = np.where(is_h)[0]
rand5 = rng.choice(hcells, size=5, replace=False)
roots = list(top10) + list(rand5)
from scipy.stats import spearmanr
rows=[]
for k, r in enumerate(roots):
    mg.uns['iroot'] = int(r)
    sc.tl.dpt(mg, n_dcs=12)
    v = mg.obs['dpt_pseudotime'].values
    ok = np.isfinite(v)&np.isfinite(base)
    rho = spearmanr(v[ok], base[ok]).statistic
    rows.append(dict(root_idx=int(r), kind='top10_homeostatic' if k<10 else 'random_homeostatic',
                     spearman_rho=rho))
    print(k, int(r), round(rho,4), flush=True)
res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/root_sensitivity.csv', index=False)
print('median |rho|:', res.spearman_rho.abs().median())
print('min |rho|:', res.spearman_rho.abs().min())
