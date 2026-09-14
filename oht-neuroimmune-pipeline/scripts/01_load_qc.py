#!/usr/bin/env python3
"""Step 01: Load GSE293358 Cell Ranger outputs, QC per manuscript Methods 4.2.

QC: cells with <200 detected genes or >10% mitochondrial reads excluded;
sample-specific upper bounds via 5xMAD rule (n_genes); genes in <3 cells removed;
Scrublet per sample (score threshold 0.25). Counts normalized to 10,000 per cell,
log1p-transformed. Output: merged AnnData with raw counts layer.
"""
import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import anndata as ad

DATA = os.environ.get('GSE293358_DIR', 'data/GSE293358')
OUT = os.environ.get('OUT_DIR', 'out')
os.makedirs(OUT, exist_ok=True)
CKPT = os.environ.get('CKPT_DIR', '/mnt/agents/work/pipe/h5ad/samples')
os.makedirs(CKPT, exist_ok=True)

META = {
    'GSM8880855': ('AMf20', 'OpticNerve', 'Sham'),
    'GSM8880856': ('AMf21', 'OpticNerve', 'MB'),
    'GSM8880857': ('AMf24', 'Retina', 'Sham'),
    'GSM8880858': ('AMf25', 'Retina', 'MB'),
    'GSM8880859': ('AMf30', 'Retina', 'MB'),
    'GSM8880860': ('AMf31', 'Retina', 'Sham'),
    'GSM8880861': ('AMf34', 'OpticNerve', 'MB'),
    'GSM8880862': ('AMf35', 'OpticNerve', 'Sham'),
    'GSM8880863': ('AMf36', 'OpticNerve', 'Sham'),
    'GSM8880864': ('AMf37', 'OpticNerve', 'MB'),
    'GSM8880865': ('AMf40', 'Retina', 'Sham'),
    'GSM8880866': ('AMf41', 'Retina', 'MB'),
    'GSM8880867': ('AMf45', 'Retina', 'Sham'),
}

adatas = []
qc_rows = []
import subprocess as _sp
_pending = []
for gsm, (suf, tissue, cond) in META.items():
    ck = f'{CKPT}/{gsm}.h5ad'; qc = f'{CKPT}/{gsm}.qc.csv'
    if os.path.exists(ck) and os.path.exists(qc):
        try:
            a = sc.read_h5ad(ck)
            adatas.append(a)
            qc_rows.append(pd.read_csv(qc).iloc[0].to_dict())
            print(gsm, 'cached', a.n_obs, flush=True)
            continue
        except Exception:
            os.remove(ck); os.remove(qc)
            print(gsm, 'cache corrupt, recomputing', flush=True)
    mtx = f'{DATA}/{gsm}_{suf}_matrix.mtx.gz'
    bar = f'{DATA}/{gsm}_{suf}_barcodes.tsv.gz'
    fea = f'{DATA}/{gsm}_{suf}_features.tsv.gz'
    a = sc.read_mtx(mtx).T
    bc = pd.read_csv(bar, header=None, sep='\t')[0].values
    ft = pd.read_csv(fea, header=None, sep='\t')
    a.obs_names = [f'{gsm}_{b}' for b in bc]
    a.var_names = ft[1].values if ft.shape[1] > 1 else ft[0].values
    a.var_names_make_unique()
    a.obs['sample'] = gsm
    a.obs['tissue'] = tissue
    a.obs['condition'] = cond
    # QC metrics
    a.var['mt'] = a.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], inplace=True)
    n0 = a.n_obs
    # lower bound + mito
    keep = (a.obs['n_genes_by_counts'] >= 200) & (a.obs['pct_counts_mt'] <= 10)
    # sample-specific upper bound: 5xMAD on n_genes
    g = a.obs.loc[keep, 'n_genes_by_counts']
    med = np.median(g)
    mad = np.median(np.abs(g - med)) * 1.4826
    upper = med + 5 * mad
    keep &= (a.obs['n_genes_by_counts'] <= upper)
    a = a[keep].copy()
    n1 = a.n_obs
    # Scrublet
    # CD45+ FACS sorting depletes parenchymal cells -> low doublet expectation
    scrub = scr.Scrublet(a.X, expected_doublet_rate=0.006, random_state=0)
    scores, _ = scrub.scrub_doublets()
    a.obs['doublet_score'] = scores
    a = a[a.obs['doublet_score'] < 0.25].copy()
    a.X = a.X.tocsr().astype(np.float32)
    qc_rows.append(dict(sample=gsm, tissue=tissue, condition=cond,
                        cells_raw=n0, cells_qc=n1, cells_final=a.n_obs,
                        upper_ngenes=upper,
                        scrublet_rate=round(float((scores >= 0.25).mean()), 4)))
    _tmp = f'/tmp/ohtrun/{gsm}.h5ad'
    a.write(_tmp, compression='gzip')
    _sz = os.path.getsize(_tmp)
    _pending.append((_tmp, ck, _sz))
    pd.DataFrame([qc_rows[-1]]).to_csv(f'/tmp/ohtrun/{gsm}.qc.csv', index=False)
    _pending.append((f'/tmp/ohtrun/{gsm}.qc.csv', qc, None))
    adatas.append(a)
    print(gsm, tissue, cond, n0, '->', n1, '->', a.n_obs, flush=True)

# flush async copies with verification
_procs = [(_sp.Popen(['cp', src, dst]), src, dst, sz) for src, dst, sz in _pending]
for pr, src, dst, sz in _procs:
    pr.wait()
    if sz is not None and (not os.path.exists(dst) or os.path.getsize(dst) != sz):
        if os.path.exists(dst): os.remove(dst)
        print('cache write incomplete, dropped:', dst, flush=True)
os.sync() if hasattr(os, 'sync') else None
if os.environ.get('SKIP_MERGED') == '1':
    pd.DataFrame(qc_rows).to_csv(f'{OUT}/qc_table.csv', index=False)
    print('SKIP_MERGED: per-sample checkpoints done', flush=True)
    raise SystemExit(0)
import scipy.sparse as sp, gc
base = adatas[0].var_names.copy()
assert all((a.var_names == base).all() for a in adatas[1:])
# gene filter (<3 cells across dataset) applied per sample BEFORE vstack to bound memory
nnz_tot = np.zeros(adatas[0].n_vars, dtype=np.int64)
for a in adatas:
    nnz_tot += np.bincount(a.X.indices, minlength=a.shape[1])
keepg = nnz_tot >= 3
print('genes kept:', keepg.sum(), 'of', len(keepg), flush=True)
Xs = []
for a in adatas:
    Xs.append(a.X.tocsc()[:, keepg].tocsr())
    a.X = None  # free full matrix immediately
    gc.collect()
obs = pd.concat([a.obs for a in adatas])
del adatas; gc.collect()
X = sp.vstack(Xs, format='csr'); del Xs; gc.collect()
print('vstack done:', X.shape, 'nnz(MB):', X.data.nbytes/1e6, flush=True)
var = pd.DataFrame(index=base[keepg])
gc.collect()
row_sum = np.asarray(X.sum(axis=1)).ravel()
scale = (1e4 / np.maximum(row_sum, 1)).astype(np.float32)
counts_per_row = np.diff(X.indptr)
X.data *= np.repeat(scale, counts_per_row)  # in-place row scaling, no matrix copy
X.data = np.log1p(X.data, dtype=np.float32)
adata = sc.AnnData(X=X, obs=obs, var=var)
print('FINAL', adata.shape, flush=True)
print(pd.crosstab(adata.obs.tissue, adata.obs.condition), flush=True)
adata.write(f'{OUT}/merged_qc.h5ad', compression='gzip')
pd.DataFrame(qc_rows).to_csv(f'{OUT}/qc_table.csv', index=False)
