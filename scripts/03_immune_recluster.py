#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
阶段2-pass2：去除非免疫污染簇 -> 免疫细胞重聚类 -> 签名打分
"""
import os, gc, json, resource
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/tmp/results/02_cluster"
os.makedirs(OUT, exist_ok=True)
sc.settings.figdir = OUT

def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
def log(m): print(m, flush=True)

adata = sc.read_h5ad("/tmp/GSE293358_clustered.h5ad")
BAD = ["14", "16", "17", "23"]   # 神经元/感光细胞污染
adata = adata[~adata.obs["leiden"].isin(BAD)].copy()
log(f"去污染后: {adata.shape} RSS {rss():.2f}")
del adata.uns["leiden"], adata.obsm["X_pca"], adata.obsm["X_umap"]
if "leiden_colors" in adata.uns: del adata.uns["leiden_colors"]
gc.collect()

# 重选 HVG -> 稀疏 PCA -> UMAP -> Leiden
sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat",
                            batch_key="sample", subset=False)
hvg_mask = adata.var["highly_variable"].to_numpy()
blocks = [adata.X[i:i + 4000][:, hvg_mask] for i in range(0, adata.n_obs, 4000)]
X_hvg = sp.vstack(blocks, format="csr").astype(np.float32)
del blocks; gc.collect()
ah = sc.AnnData(X=X_hvg, obs=adata.obs.copy(), var=adata.var[hvg_mask].copy())
sc.pp.scale(ah, zero_center=False, max_value=10)
sc.pp.pca(ah, n_comps=50, zero_center=False, svd_solver="arpack")
adata.obsm["X_pca"] = ah.obsm["X_pca"].astype(np.float32)
del ah, X_hvg; gc.collect()
log(f"PCA 完成 RSS {rss():.2f}")

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden2")
log(f"重聚类: {adata.obs['leiden2'].nunique()} 簇 RSS {rss():.2f}")

MARKERS = {
    "Microglia":     ["P2ry12", "Tmem119", "Sall1", "Hexb", "Siglech", "Cx3cr1"],
    "DAM_Apoe":      ["Apoe", "Ctsb", "Ctsd", "Lyz2", "Cd63", "C1qb"],
    "MHCII":         ["Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1"],
    "Mono/Macro":    ["Ccr2", "Ly6c2", "Plac8", "Chil3", "Ms4a7"],
    "BAM":           ["Mrc1", "Lyve1", "Pf4", "F13a1"],
    "DC":            ["Flt3", "Itgax", "Xcr1", "Clec9a", "Cd209a"],
    "T":             ["Cd3e", "Cd3g", "Cd8a", "Cd4", "Trac", "Il7r"],
    "NK":            ["Ncr1", "Klrk1", "Gzma", "Nkg7"],
    "B":             ["Cd79a", "Cd79b", "Ms4a1", "Pax5"],
    "Neutrophil":    ["S100a8", "S100a9", "Ly6g", "Cxcl2", "Il1b"],
    "Prolif":        ["Mki67", "Top2a", "Stmn1"],
}
rows = []
for ct, genes in MARKERS.items():
    g_use = [g for g in genes if g in adata.var_names]
    sub = adata[:, g_use].X
    if sp.issparse(sub): sub = sub.toarray()
    df = pd.DataFrame(sub, columns=g_use, index=adata.obs_names)
    df["leiden2"] = adata.obs["leiden2"].to_numpy()
    m = df.groupby("leiden2").mean().mean(axis=1)
    for cl, v in m.items():
        rows.append({"cluster": cl, "sig": ct, "expr": v})
piv = pd.DataFrame(rows).pivot_table(index="cluster", columns="sig", values="expr")
piv.to_csv(os.path.join(OUT, "cluster2_signature_scores.csv"))
print(piv.round(2).to_string(), flush=True)

comp = pd.crosstab(adata.obs["leiden2"], [adata.obs["tissue"], adata.obs["condition"]])
comp.to_csv(os.path.join(OUT, "cluster2_composition_counts.csv"))
print(comp.to_string(), flush=True)

sc.pl.umap(adata, color=["leiden2"], legend_loc="on data", show=False, frameon=False)
plt.savefig(os.path.join(OUT, "UMAP_leiden2.png"), dpi=200, bbox_inches="tight"); plt.close()

adata.write_h5ad("/tmp/GSE293358_immune.h5ad", compression="gzip")
log(f"保存完成 RSS {rss():.2f}")
print(json.dumps({"cells": int(adata.n_obs), "clusters": int(adata.obs["leiden2"].nunique())}), flush=True)
