#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""阶段3：小胶质亚群重聚类 + 扩散拟时序（DPT）+ 轨迹基因动态"""
import os, gc, json, resource
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/tmp/results/03_mg"
os.makedirs(OUT, exist_ok=True)
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
def log(m): print(m, flush=True)

adata = sc.read_h5ad("/tmp/GSE293358_final.h5ad")
MG_TYPES = ["Homeostatic MG", "MHCII+ MG", "Activated MG (Ccl5+)",
            "DAM-like MG", "Iron+ MG (Fth1+)", "Proliferating MG"]
mg = adata[adata.obs["celltype"].isin(MG_TYPES)].copy()
del adata; gc.collect()
log(f"小胶质细胞: {mg.shape} RSS {rss():.2f}")

# ---------- 1. 小胶质内重聚类 ----------
for k in ["leiden2", "X_pca", "X_umap"]:
    if k in mg.obsm: del mg.obsm[k]
sc.pp.highly_variable_genes(mg, n_top_genes=2000, flavor="seurat",
                            batch_key="sample", subset=False)
hvg = mg.var["highly_variable"].to_numpy()
blocks = [mg.X[i:i+4000][:, hvg] for i in range(0, mg.n_obs, 4000)]
Xh = sp.vstack(blocks, format="csr").astype(np.float32)
del blocks; gc.collect()
ah = sc.AnnData(X=Xh, obs=mg.obs.copy(), var=mg.var[hvg].copy())
sc.pp.scale(ah, zero_center=False, max_value=10)
sc.pp.pca(ah, n_comps=40, zero_center=False, svd_solver="arpack")
mg.obsm["X_pca"] = ah.obsm["X_pca"].astype(np.float32)
del ah, Xh; gc.collect()
sc.pp.neighbors(mg, n_neighbors=15, n_pcs=25)
sc.tl.umap(mg)
sc.tl.leiden(mg, resolution=0.8, key_added="mg_state")
log(f"MG 亚群: {mg.obs['mg_state'].nunique()} RSS {rss():.2f}")

# ---------- 2. 状态签名打分 ----------
SIG = {
    "Homeostatic": ["P2ry12", "Tmem119", "Sall1", "Selplg", "Siglech"],
    "DAM":         ["Apoe", "Ctsb", "Ctsd", "Cd63", "Lyz2", "Trem2", "Tyrobp", "Axl", "Cst7", "Spp1", "Lpl"],
    "MHCII":       ["Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1"],
    "IFN":         ["Ifit1", "Ifit3", "Isg15", "Irf7", "Stat1", "Ifi204"],
    "Chemokine":   ["Ccl5", "Ccl4", "Ccl12", "Cxcl10", "Cxcl16"],
    "Iron":        ["Fth1", "Ftl1"],
    "Prolif":      ["Mki67", "Top2a", "Stmn1"],
    "Ribo":        ["Rps12", "Rps24", "Rpl13", "Rpl37a"],
}
rows = []
for sname, genes in SIG.items():
    g_use = [g for g in genes if g in mg.var_names]
    sub = mg[:, g_use].X
    if sp.issparse(sub): sub = sub.toarray()
    df = pd.DataFrame(sub.mean(axis=1), columns=[sname])
    df["mg_state"] = mg.obs["mg_state"].to_numpy()
    rows.append(df.groupby("mg_state")[sname].mean())
sig_tab = pd.concat(rows, axis=1)
sig_tab.to_csv(os.path.join(OUT, "mg_state_signatures.csv"))
print(sig_tab.round(2).to_string(), flush=True)

# 组成
comp = pd.crosstab(mg.obs["mg_state"], [mg.obs["tissue"], mg.obs["condition"]])
comp.to_csv(os.path.join(OUT, "mg_state_composition.csv"))
print(comp.to_string(), flush=True)

# ---------- 3. 扩散图 + DPT 拟时序 ----------
sc.tl.diffmap(mg, n_comps=12)
# 用稳态得分最高的细胞做根
hgenes = [g for g in ["P2ry12", "Tmem119", "Sall1"] if g in mg.var_names]
hscore = np.asarray(mg[:, hgenes].X.mean(axis=1)).ravel()
mg.uns["iroot"] = int(np.argmax(hscore))
sc.tl.dpt(mg, n_dcs=12)
log(f"DPT 完成 RSS {rss():.2f}")

# ---------- 4. 轨迹相关基因（对 HVG 做 spearman 相关） ----------
dpt = mg.obs["dpt_pseudotime"].to_numpy()
valid = np.isfinite(dpt)
Xh = mg[:, hvg].X
if sp.issparse(Xh): Xh = Xh.toarray()
Xv = Xh[valid]; dv = dpt[valid]
# 抽样至多 8000 细胞加速
if len(dv) > 8000:
    idx = np.random.RandomState(0).choice(len(dv), 8000, replace=False)
    Xv = Xv[idx]; dv = dv[idx]
res = []
for j in range(Xv.shape[1]):
    rho, p = stats.spearmanr(dv, Xv[:, j])
    res.append((mg.var_names[hvg][j], rho))
traj = pd.DataFrame(res, columns=["gene", "spearman_dpt"]).sort_values("spearman_dpt")
traj.to_csv(os.path.join(OUT, "traj_gene_correlations.csv"), index=False)
log("轨迹相关基因完成")
print("负相关 top15（激活时丢失）:", ", ".join(traj.head(15)["gene"]), flush=True)
print("正相关 top15（激活时获得）:", ", ".join(traj.tail(15)["gene"]), flush=True)

mg.write_h5ad("/tmp/GSE293358_mg.h5ad", compression="gzip")
log("保存 mg.h5ad")
print(json.dumps({"mg_cells": int(mg.n_obs), "states": int(mg.obs['mg_state'].nunique())}), flush=True)
