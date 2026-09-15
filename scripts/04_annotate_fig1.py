#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""阶段2-pass3：终注释 + Fig1 图谱出图"""
import os, gc, json, resource
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/tmp/results/02_cluster"
os.makedirs(OUT, exist_ok=True)
sc.settings.figdir = OUT
def log(m): print(m, flush=True)

LABEL = {
    "0": "Homeostatic MG", "1": "Homeostatic MG", "2": "Homeostatic MG",
    "3": "Homeostatic MG", "4": "Homeostatic MG", "5": "Homeostatic MG",
    "6": "Homeostatic MG", "8": "Homeostatic MG",
    "7": "MHCII+ MG",
    "9": "Activated MG (Ccl5+)",
    "10": "DAM-like MG", "12": "DAM-like MG",
    "18": "Iron+ MG (Fth1+)",
    "17": "Proliferating MG",
    "14": "BAM",
    "15": "MHCII+ Mac/APC",
    "19": "pDC",
    "13": "T cell", "21": "gdT cell",
    "16": "NK cell", "20": "B cell",
    "11": "Neutrophil",
}

adata = sc.read_h5ad("/tmp/GSE293358_immune.h5ad")
adata = adata[adata.obs["leiden2"] != "22"].copy()   # 剔除 13 细胞杂簇
adata.obs["celltype"] = adata.obs["leiden2"].map(LABEL).astype("category")
log(f"注释完成: {adata.shape}, {adata.obs['celltype'].nunique()} 类")
print(adata.obs["celltype"].value_counts().to_string(), flush=True)

# ---- Fig1a: UMAP by celltype ----
ct_order = ["Homeostatic MG", "MHCII+ MG", "Activated MG (Ccl5+)", "DAM-like MG",
            "Iron+ MG (Fth1+)", "Proliferating MG", "BAM", "MHCII+ Mac/APC",
            "pDC", "T cell", "gdT cell", "NK cell", "B cell", "Neutrophil"]
adata.obs["celltype"] = adata.obs["celltype"].cat.set_categories(ct_order)
palette = ["#205B57", "#3E7A6E", "#5C988B", "#8FA8A0", "#B5C8B8", "#D4A574",
           "#6B3A5B", "#8E6278", "#C9A0B5", "#37455C", "#546E7A", "#B8860B",
           "#7C4D3E", "#C0532F"]
sc.pl.umap(adata, color="celltype", palette=palette, show=False, frameon=False,
           size=3, legend_fontsize=8)
plt.savefig(os.path.join(OUT, "Fig1a_UMAP_celltype.png"), dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig1b: UMAP split by condition ----
sc.pl.umap(adata, color="condition", groups=["Sham", "MB"], show=False, frameon=False,
           size=3, palette={"Sham": "#8FA8A0", "MB": "#C0532F"}, ncols=2)
plt.savefig(os.path.join(OUT, "Fig1b_UMAP_condition.png"), dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig1c: marker dotplot ----
MK = {
    "Homeostatic MG": ["P2ry12", "Tmem119", "Sall1"],
    "MHCII+ MG": ["H2-Aa", "H2-Ab1"],
    "Activated MG (Ccl5+)": ["Ccl5", "Apoe"],
    "DAM-like MG": ["Cd63", "Ctsb", "Ctsd"],
    "Iron+ MG (Fth1+)": ["Fth1", "Ftl1"],
    "Proliferating MG": ["Mki67", "Top2a"],
    "BAM": ["Mrc1", "Ms4a7"],
    "MHCII+ Mac/APC": ["Cd74", "H2-Eb1"],
    "pDC": ["Tcf4", "Bst2"],
    "T cell": ["Cd3e", "Cd8a"],
    "gdT cell": ["Trdc", "Il7r"],
    "NK cell": ["Nkg7", "Ncr1"],
    "B cell": ["Cd79a", "Ms4a1"],
    "Neutrophil": ["S100a8", "S100a9"],
}
mk_flat = [g for v in MK.values() for g in v if g in adata.var_names]
sc.pl.dotplot(adata, mk_flat, groupby="celltype", show=False, dendrogram=False,
              categories_order=ct_order, cmap="viridis")
plt.savefig(os.path.join(OUT, "Fig1c_dotplot.png"), dpi=300, bbox_inches="tight"); plt.close()
log("Fig1a-c 完成")

# ---- Fig1d: 组成堆叠条 + 每重复比例统计 ----
tab = pd.crosstab([adata.obs["sample"], adata.obs["tissue"], adata.obs["condition"]],
                  adata.obs["celltype"])
prop = tab.div(tab.sum(axis=1), axis=0)
prop.to_csv(os.path.join(OUT, "celltype_proportions_per_sample.csv"))

grp = pd.crosstab([adata.obs["tissue"], adata.obs["condition"]], adata.obs["celltype"])
grp_prop = grp.div(grp.sum(axis=1), axis=0).T.loc[ct_order]
fig, ax = plt.subplots(figsize=(7, 5))
bottom = np.zeros(grp_prop.shape[1])
for ct, color in zip(grp_prop.index, palette):
    ax.bar(range(grp_prop.shape[1]), grp_prop.loc[ct], bottom=bottom,
           label=ct, color=color, width=0.65, edgecolor="white", linewidth=0.3)
    bottom += grp_prop.loc[ct].to_numpy()
ax.set_xticks(range(grp_prop.shape[1]))
ax.set_xticklabels([f"{t}\n{c}" for t, c in grp_prop.columns], fontsize=9)
ax.set_ylabel("Proportion"); ax.set_ylim(0, 1)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Fig1d_composition.png"), dpi=300, bbox_inches="tight"); plt.close()

# ---- 统计: 每个组织内 MB vs Sham 各细胞类型比例 Welch t 检验 ----
rows = []
prop_df = prop.reset_index()
for tissue in ["Retina", "OpticNerve"]:
    sub = prop_df[prop_df["tissue"] == tissue]
    for ct in ct_order:
        a = sub.loc[sub["condition"] == "MB", ct].to_numpy()
        b = sub.loc[sub["condition"] == "Sham", ct].to_numpy()
        if a.sum() == 0 and b.sum() == 0: continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        lfc = np.log2((a.mean() + 1e-4) / (b.mean() + 1e-4))
        rows.append({"tissue": tissue, "celltype": ct,
                     "prop_MB": a.mean(), "prop_Sham": b.mean(),
                     "log2FC": lfc, "pvalue": p})
st = pd.DataFrame(rows)
st["padj"] = st.groupby("tissue")["pvalue"].transform(
    lambda p: np.minimum(p * len(p) / (np.argsort(np.argsort(p)) + 1), 1))
st.to_csv(os.path.join(OUT, "celltype_proportion_stats.csv"), index=False)
print(st.round(4).to_string(), flush=True)

adata.write_h5ad("/tmp/GSE293358_final.h5ad", compression="gzip")
log("保存最终 h5ad 完成")
print(json.dumps({"cells": int(adata.n_obs), "celltypes": int(adata.obs['celltype'].nunique())}), flush=True)
