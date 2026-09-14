#!/usr/bin/env python3
"""Step 03: LIANA cell-cell communication, one tissue x condition subset per process."""
import sys, os, numpy as np, pandas as pd, scanpy as sc, liana as li

DATA = os.environ.get('GSE293358_DIR', 'data/GSE293358')
OUT = os.environ.get('OUT_DIR', 'out')
tissue, cond, nperms = sys.argv[1], sys.argv[2], int(sys.argv[3])

adata = sc.read_h5ad(f'{OUT}/annotated.h5ad')
mask = (adata.obs['tissue'] == tissue) & (adata.obs['condition'] == cond) & (~adata.obs['is_contaminant'])
sub = adata[mask].copy()
del adata
print(f'== {tissue} {cond}: {sub.n_obs} cells, {sub.obs["cell_class"].nunique()} classes ==', flush=True)

li.method.rank_aggregate(adata=sub, groupby='cell_class', resource_name='mouseconsensus',
                         expr_prop=0.1, min_cells=10, n_perms=nperms, seed=1337,
                         use_raw=False, verbose=True)
res = sub.uns['liana_res'].copy()
res['tissue'], res['condition'], res['n_perms'] = tissue, cond, nperms
res.to_csv(f'{OUT}/liana_{tissue}_{cond}_{nperms}perms.csv.gz', index=False, compression='gzip')
print(f'{tissue} {cond}: {res.shape[0]} interactions -> saved', flush=True)

spp = res[(res.ligand_complex == 'Spp1') & (res.receptor_complex == 'Cd44')].sort_values('lr_means', ascending=False)
print(f'\nSpp1-Cd44 rows ({tissue} {cond}):', flush=True)
print(spp[['source','target','lr_means','specificity_rank','magnitude_rank']].head(12).to_string(index=False), flush=True)
