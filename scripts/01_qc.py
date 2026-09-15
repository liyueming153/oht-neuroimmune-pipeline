#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
阶段1：GSE293358 单细胞数据读入与质控（内存安全版，cgroup 上限 3GB）
策略：逐样本完成 QC 指标与细胞过滤 -> 全局基因计数 -> 逐样本基因过滤 -> 合并
"""
import os, re, gzip, glob, json, gc, resource
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/mnt/agents/output/oht_project"
RAW = "/tmp/gse_raw"
OUT = "/tmp/results/01_qc"
os.makedirs(OUT, exist_ok=True)

def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

SAMPLE_META = {
    "GSM8880855": ("OpticNerve", "Sham", "S1"),
    "GSM8880856": ("OpticNerve", "MB",   "S1"),
    "GSM8880857": ("Retina",     "Sham", "S1"),
    "GSM8880858": ("Retina",     "MB",   "S1"),
    "GSM8880859": ("Retina",     "MB",   "S2"),
    "GSM8880860": ("Retina",     "Sham", "S2"),
    "GSM8880861": ("OpticNerve", "MB",   "S2"),
    "GSM8880862": ("OpticNerve", "Sham", "S2"),
    "GSM8880863": ("OpticNerve", "Sham", "S3"),
    "GSM8880864": ("OpticNerve", "MB",   "S3"),
    "GSM8880865": ("Retina",     "Sham", "S3"),
    "GSM8880866": ("Retina",     "MB",   "S3"),
    "GSM8880867": ("Retina",     "Sham", "S4"),
}

def load_sample(gsm):
    files = glob.glob(os.path.join(RAW, gsm + "_*"))
    mtx_f = [f for f in files if f.endswith("matrix.mtx.gz")][0]
    bar_f = [f for f in files if "barcodes" in f][0]
    fea_f = [f for f in files if "features" in f or "genes" in f][0]
    X = sio.mmread(gzip.open(mtx_f, "rb")).tocsr().T.tocsr().astype(np.int32)
    with gzip.open(bar_f, "rt") as fh:
        barcodes = [l.strip().split("\t")[0] for l in fh]
    with gzip.open(fea_f, "rt") as fh:
        genes = [l.strip().split("\t") for l in fh]
    gene_names = [g[1] if len(g) > 1 else g[0] for g in genes]
    adata = sc.AnnData(X=X)
    adata.obs_names = [f"{gsm}_{b}" for b in barcodes]
    adata.var_names = gene_names
    adata.var_names_make_unique()
    return adata

def mad_cut(x, n=5):
    med = np.median(x); mad = np.median(np.abs(x - med))
    return med + n * 1.4826 * mad

# ---------- 1. 逐样本：读入 -> QC 指标 -> 细胞过滤 ----------
adatas = []
obs_pre_all = []   # 质控前指标（绘图用）
gene_nn_global = None
gene_names_ref = None
n_raw = 0
for gsm, (tissue, cond, rep) in SAMPLE_META.items():
    a = load_sample(gsm)
    if gene_names_ref is None:
        gene_names_ref = a.var_names.to_numpy()
        mt_mask = np.asarray(pd.Index(gene_names_ref).str.startswith(("mt-", "Mt-", "MT-")))
        gene_nn_global = np.zeros(len(gene_names_ref), dtype=np.int64)
    a.var["mt"] = mt_mask
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, log1p=False)
    gene_nn_global += np.asarray(a.X.getnnz(axis=0)).ravel()
    n_raw += a.n_obs

    ob = a.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    ob["sample"] = gsm; ob["tissue"] = tissue; ob["condition"] = cond
    obs_pre_all.append(ob)

    upper = mad_cut(a.obs["n_genes_by_counts"].to_numpy())
    keep = (a.obs["n_genes_by_counts"] >= 200) & \
           (a.obs["n_genes_by_counts"] <= upper) & \
           (a.obs["pct_counts_mt"] <= 10)
    a = a[keep.to_numpy()].copy()
    a.obs["sample"] = gsm; a.obs["tissue"] = tissue
    a.obs["condition"] = cond; a.obs["replicate"] = rep
    print(f"{gsm} {tissue}-{cond}-{rep}: {keep.sum()}/{len(keep)} 细胞保留 (上限 {upper:.0f}) RSS {rss():.2f}", flush=True)
    adatas.append(a)
    del ob; gc.collect()

obs_pre = pd.concat(obs_pre_all); del obs_pre_all; gc.collect()
print(f"原始总细胞 {n_raw}, 过滤后 {sum(a.n_obs for a in adatas)}", flush=True)

# ---------- 2. 全局基因过滤（min_cells=3），逐样本列切片避免大拷贝 ----------
gene_keep = gene_nn_global >= 3
print(f"基因过滤: {len(gene_keep)} -> {int(gene_keep.sum())}", flush=True)
adatas = [a[:, gene_keep].copy() for a in adatas]
gc.collect()

# ---------- 3. 质控前小提琴图（基于 obs_pre，纯 pandas 绘图） ----------
def violin3(df, path, group):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))
    metrics = [("n_genes_by_counts", "Genes"), ("total_counts", "UMI"), ("pct_counts_mt", "mt%")]
    cats = list(dict.fromkeys(df[group]))
    for ax, (m, lab) in zip(axes, metrics):
        data = [df.loc[df[group] == c, m].to_numpy() for c in cats]
        vp = ax.violinplot(data, showmedians=True)
        for b in vp["bodies"]:
            b.set_facecolor("#205B57"); b.set_alpha(0.7)
        for k in ("cmedians", "cmins", "cmaxes", "cbars"):
            if k in vp: vp[k].set_color("#3E2723")
        ax.set_xticks(range(1, len(cats) + 1)); ax.set_xticklabels(cats, rotation=90, fontsize=7)
        ax.set_title(lab)
    plt.tight_layout(); plt.savefig(path, dpi=200); plt.close()

violin3(obs_pre, os.path.join(OUT, "QC_violin_prefilter.png"), "sample")
print("prefilter 图完成 RSS", round(rss(), 2), flush=True)

# ---------- 4. 合并 ----------
adata = sc.concat(adatas, join="outer", index_unique=None)
adata.var_names_make_unique()
del adatas; gc.collect()
print("合并后:", adata.shape, "RSS", round(rss(), 2), flush=True)

# 质控后图（按 condition 分面，从 obs 画）
violin3(adata.obs, os.path.join(OUT, "QC_violin_postfilter.png"), "condition")
print("postfilter 图完成", flush=True)

# ---------- 5. 双联体检测（每样本 Scrublet） ----------
try:
    import scrublet as scr
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    for gsm in adata.obs["sample"].unique():
        mask = (adata.obs["sample"] == gsm).to_numpy()
        sub_X = adata.X[mask].tocsr()
        scrub = scr.Scrublet(sub_X)
        score, pred = scrub.scrub_doublets(verbose=False)
        idx = np.where(mask)[0]
        adata.obs.iloc[idx, adata.obs.columns.get_loc("doublet_score")] = score
        adata.obs.iloc[idx, adata.obs.columns.get_loc("predicted_doublet")] = pred
        del sub_X, scrub; gc.collect()
        print(f"  Scrublet {gsm}: {int(pred.sum())}/{len(pred)} 双联体 RSS {rss():.2f}", flush=True)
    nd = int(adata.obs["predicted_doublet"].sum())
    print(f"Scrublet 预测双联体: {nd} ({nd / adata.n_obs * 100:.1f}%)", flush=True)
    keep_mask = ~adata.obs["predicted_doublet"].to_numpy()
    adata2 = adata[keep_mask].copy()
    del adata, keep_mask; gc.collect()
    adata = adata2; del adata2
    print("去双联体后:", adata.shape, "RSS", round(rss(), 2), flush=True)
except Exception as e:
    print("[WARN] Scrublet 跳过:", repr(e), flush=True)

# ---------- 6. 汇总与保存 ----------
summary = (adata.obs.groupby(["sample", "tissue", "condition", "replicate"])
           .size().reset_index(name="cells_postQC"))
summary.to_csv(os.path.join(OUT, "cells_per_sample_postQC.csv"), index=False)

fig, ax = plt.subplots(figsize=(8, 4))
piv = summary.pivot_table(index="sample", columns="condition", values="cells_postQC")
piv.plot(kind="bar", ax=ax, color=["#8FA8A0", "#205B57"])
ax.set_ylabel("Cells post-QC"); plt.xticks(rotation=45, ha="right")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "cells_per_sample.png"), dpi=200); plt.close()

adata.write_h5ad("/tmp/GSE293358_qc.h5ad", compression="gzip")
print("QC 完成，输出目录:", OUT, flush=True)
print(json.dumps({"raw_cells": int(n_raw), "final_cells": int(adata.n_obs),
                  "final_genes": int(adata.n_vars)}, ensure_ascii=False), flush=True)
