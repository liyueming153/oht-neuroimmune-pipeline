#!/usr/bin/env python3
"""Step 02a: dimensionality reduction and Leiden clustering (Methods 4.2).

Input is already normalized (log1p of 1e4-normalized counts) by step 01.
Embedding built on 2500 HVGs (seurat, batch_key=sample):
scale(zero_center=False) -> PCA(50) -> kNN(k=15, 30 PCs) -> UMAP ->
Leiden(res=1.0, key 'leiden1'). Labels and UMAP are transferred back to
the full-gene object so downstream marker scoring uses all genes.
"""
import os
import numpy as np
import scanpy as sc

OUT = os.environ.get('OUT_DIR', 'out')
CKPT = os.environ.get('CKPT_DIR', '/mnt/agents/work/pipe/h5ad/samples')
if os.path.exists(f'{OUT}/merged_qc.h5ad'):
    full = sc.read_h5ad(f'{OUT}/merged_qc.h5ad')
else:  # assemble from per-sample checkpoints (identical to 01 tail)
    import glob, scipy.sparse as sp, pandas as pd, gc
    files = sorted(glob.glob(f'{CKPT}/*.h5ad'))
    assert len(files) == 13, files
    ads = [sc.read_h5ad(f) for f in files]
    base = ads[0].var_names.copy()
    nnz_tot = np.zeros(ads[0].n_vars, dtype=np.int64)
    for a in ads: nnz_tot += np.bincount(a.X.indices, minlength=a.shape[1])
    keepg = nnz_tot >= 3
    Xs = [a.X.tocsc()[:, keepg].tocsr() for a in ads]
    obs = pd.concat([a.obs for a in ads])
    X = sp.vstack(Xs, format='csr'); del Xs, ads; gc.collect()
    row_sum = np.asarray(X.sum(axis=1)).ravel()
    scale = (1e4 / np.maximum(row_sum, 1)).astype(np.float32)
    X.data *= np.repeat(scale, np.diff(X.indptr))
    X.data = np.log1p(X.data, dtype=np.float32)
    full = sc.AnnData(X=X, obs=obs, var=pd.DataFrame(index=base[keepg]))
    del X, obs; gc.collect()
print('loaded', full.shape, flush=True)

sc.pp.highly_variable_genes(full, n_top_genes=2500, flavor='seurat', batch_key='sample', subset=False)
emb = full[:, full.var.highly_variable].copy()
sc.pp.scale(emb, zero_center=False)
sc.pp.pca(emb, n_comps=50)
sc.pp.neighbors(emb, n_neighbors=15, n_pcs=30)
sc.tl.umap(emb)
sc.tl.leiden(emb, resolution=1.0, key_added='leiden1', flavor='igraph', n_iterations=2)
print(emb.obs['leiden1'].value_counts().sort_index(), flush=True)

full.obs['leiden1'] = emb.obs['leiden1']
full.obsm['X_umap'] = emb.obsm['X_umap']
full.write(f'{OUT}/clustered_ckpt.h5ad', compression='gzip')
print('saved clustered_ckpt.h5ad', flush=True)
from keep_genes import EXTRA_GENES
keep = full.var.highly_variable | full.var_names.isin(EXTRA_GENES)
slim = full[:, keep].copy()
slim.write(f'{OUT}/clustered_slim.h5ad', compression='gzip')
print('saved clustered_slim.h5ad', slim.shape, flush=True)
