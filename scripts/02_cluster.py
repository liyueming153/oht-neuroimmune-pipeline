#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
阶段2-pass1：标准化 -> HVG -> PCA -> UMAP -> Leiden 聚类 -> marker 打分
内存安全：cgroup 3GB，逐步骤控制峰值
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

adata = sc.read_h5ad("/tmp/GSE293358_qc.h5ad")
log(f"载入 {adata.shape} RSS {rss():.2f}")

# ---------- 1. 标准化 ----------
# 注意：不保留 counts layer（/tmp/GSE293358_qc.h5ad 里就是 int32 原始计数，阶段3可取）
adata.X = adata.X.astype(np.float32)
gc.collect()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
log(f"标准化完成 RSS {rss():.2f}")

# ---------- 2. HVG（按样本去批次） ----------
sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat",
                            batch_key="sample", subset=False)
hvg_n = int(adata.var["highly_variable"].sum())
log(f"HVG: {hvg_n} RSS {rss():.2f}")

# ---------- 3. 仅对 HVG scale + PCA（分块列切片防 CSC 转换 OOM） ----------
hvg_mask = adata.var["highly_variable"].to_numpy()
blocks = []
for i in range(0, adata.n_obs, 4000):
    blocks.append(adata.X[i:i + 4000][:, hvg_mask])
X_hvg = sp.vstack(blocks, format="csr").astype(np.float32)
del blocks; gc.collect()
adata_hvg = sc.AnnData(X=X_hvg, obs=adata.obs.copy(),
                       var=adata.var[hvg_mask].copy())
del X_hvg; gc.collect()
sc.pp.scale(adata_hvg, zero_center=False, max_value=10)  # 保持稀疏，防致密化 OOM
sc.pp.pca(adata_hvg, n_comps=50, zero_center=False, svd_solver="arpack")
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].astype(np.float32)
adata.uns["pca_variance_ratio"] = adata_hvg.uns["pca"]["variance_ratio"]
del adata_hvg; gc.collect()
log(f"PCA 完成 RSS {rss():.2f}")

# ---------- 4. 邻接图 + UMAP + Leiden ----------
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden")
n_cl = adata.obs["leiden"].nunique()
log(f"UMAP+Leiden 完成: {n_cl} 个簇 RSS {rss():.2f}")

# ---------- 5. 每簇 marker（wilcoxon） ----------
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=30)
mk = sc.get.rank_genes_groups_df(adata, group=None)
mk.to_csv(os.path.join(OUT, "cluster_markers_top30.csv"), index=False)
log("marker 计算完成")

# ---------- 6. 经典 marker 打分（细胞类型注释依据） ----------
MARKERS = {
    "Microglia":        ["P2ry12", "Tmem119", "Sall1", "Hexb", "Siglech", "Cx3cr1"],
    "Mono/Macro":       ["Ccr2", "Ly6c2", "Plac8", "Chil3", "Ms4a7"],
    "MHCII+ AP":        ["Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1"],
    "DC":               ["Flt3", "Itgax", "Xcr1", "Clec9a", "Cd209a"],
    "BAM/Perivasc":     ["Mrc1", "Lyve1", "Pf4", "F13a1"],
    "T cell":           ["Cd3e", "Cd3g", "Cd8a", "Cd4", "Trac"],
    "NK cell":          ["Ncr1", "Klrk1", "Gzma", "Nkg7"],
    "B cell":           ["Cd79a", "Cd79b", "Ms4a1", "Cd19"],
    "Neutrophil":       ["S100a8", "S100a9", "Ly6g", "Mpo"],
    "Proliferating":    ["Mki67", "Top2a", "Stmn1"],
    "RBC contam":       ["Hba-a1", "Hba-a2", "Hbb-bs"],
}
rows = []
for ct, genes in MARKERS.items():
    g_use = [g for g in genes if g in adata.var_names]
    sub = adata[:, g_use].X
    if sp.issparse(sub): sub = sub.toarray()
    df = pd.DataFrame(sub, columns=g_use, index=adata.obs_names)
    df["leiden"] = adata.obs["leiden"].to_numpy()
    m = df.groupby("leiden").mean()
    m["score"] = m.mean(axis=1)
    for cl in m.index:
        rows.append({"cluster": cl, "celltype_sig": ct, "mean_expr": m.loc[cl, "score"]})
score_df = pd.DataFrame(rows)
piv = score_df.pivot_table(index="cluster", columns="celltype_sig", values="mean_expr")
piv.to_csv(os.path.join(OUT, "cluster_signature_scores.csv"))
log("signature 打分完成")
print(piv.round(2).to_string(), flush=True)

# ---------- 7. 图 ----------
sc.pl.umap(adata, color=["leiden"], legend_loc="on data", show=False,
           title="Leiden clusters", frameon=False)
plt.savefig(os.path.join(OUT, "UMAP_leiden.png"), dpi=200, bbox_inches="tight"); plt.close()

sc.pl.umap(adata, color=["condition"], show=False, frameon=False,
           palette={"Sham": "#8FA8A0", "MB": "#C0532F"})
plt.savefig(os.path.join(OUT, "UMAP_condition.png"), dpi=200, bbox_inches="tight"); plt.close()

sc.pl.umap(adata, color=["tissue"], show=False, frameon=False,
           palette={"Retina": "#205B57", "OpticNerve": "#D4A574"})
plt.savefig(os.path.join(OUT, "UMAP_tissue.png"), dpi=200, bbox_inches="tight"); plt.close()

mk_flat = sorted({g for v in MARKERS.values() for g in v if g in adata.var_names})
sc.pl.dotplot(adata, mk_flat, groupby="leiden", show=False, dendrogram=False)
plt.savefig(os.path.join(OUT, "dotplot_markers.png"), dpi=200, bbox_inches="tight"); plt.close()
log("图完成")

# ---------- 8. 簇 x 条件/组织 组成表 ----------
comp = pd.crosstab(adata.obs["leiden"], [adata.obs["tissue"], adata.obs["condition"]])
comp.to_csv(os.path.join(OUT, "cluster_composition_counts.csv"))
print(comp.to_string(), flush=True)

adata.write_h5ad("/tmp/GSE293358_clustered.h5ad", compression="gzip")
log(f"保存完成 RSS {rss():.2f}")
print(json.dumps({"clusters": n_cl, "cells": int(adata.n_obs)}), flush=True)
