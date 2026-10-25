#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""阶段3-出图：Fig2 小胶质亚群 + Fig3 拟时序"""
import os, gc, json
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
def log(m): print(m, flush=True)

mg = sc.read_h5ad("/tmp/GSE293358_mg.h5ad")

STATE_LABEL = {
    "0": "MG-Homeo-ON", "1": "MG-Homeo-Ret", "2": "MG-Homeo-Ret2",
    "3": "MG-Homeo-ON2", "4": "MG-MHCII+", "5": "MG-lowRNA",
    "6": "MG-IFN", "7": "MG-Ccl5+ act.", "8": "MG-Prolif",
    "9": "MG-DAM", "10": "MG-MHCII+ (Ret)", "11": "MG-DAM-Iron",
}
mg.obs["mg_label"] = mg.obs["mg_state"].map(STATE_LABEL).astype("category")
order = ["MG-Homeo-ON", "MG-Homeo-Ret", "MG-Homeo-Ret2", "MG-Homeo-ON2",
         "MG-lowRNA", "MG-IFN", "MG-MHCII+", "MG-MHCII+ (Ret)",
         "MG-Ccl5+ act.", "MG-DAM", "MG-DAM-Iron", "MG-Prolif"]
mg.obs["mg_label"] = mg.obs["mg_label"].cat.set_categories(order)
pal = ["#1B4B47", "#205B57", "#3E7A6E", "#5C988B", "#B5C8B8", "#546E7A",
       "#B8860B", "#D4A574", "#8E6278", "#C0532F", "#7C4D3E", "#6B3A5B"]
mg.obs["group"] = (mg.obs["tissue"].astype(str) + "-" + mg.obs["condition"].astype(str)).astype("category")

# ---- Fig2a UMAP ----
sc.pl.umap(mg, color="mg_label", palette=pal, show=False, frameon=False, size=4,
           legend_fontsize=7)
plt.savefig(f"{OUT}/Fig2a_UMAP_mg_states.png", dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig2b UMAP by group (4 分面) ----
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
um = mg.obs[["group"]].copy()
xy = pd.DataFrame(mg.obsm["X_umap"], columns=["u1", "u2"], index=mg.obs_names)
for ax, gname, color in zip(axes.ravel(),
        ["Retina-Sham", "Retina-MB", "OpticNerve-Sham", "OpticNerve-MB"],
        ["#8FA8A0", "#C0532F", "#8FA8A0", "#C0532F"]):
    m = (mg.obs["group"] == gname).to_numpy()
    ax.scatter(xy["u1"], xy["u2"], s=0.3, c="#DDDDDD", rasterized=True)
    ax.scatter(xy.loc[m, "u1"], xy.loc[m, "u2"], s=0.5, c=color, rasterized=True)
    ax.set_title(f"{gname} (n={m.sum()})", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig2b_UMAP_by_group.png", dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig2c 组成堆叠条 ----
grp = pd.crosstab(mg.obs["group"], mg.obs["mg_label"])
grp_prop = grp.div(grp.sum(axis=1), axis=0).T.loc[order]
fig, ax = plt.subplots(figsize=(6.5, 5))
bottom = np.zeros(grp_prop.shape[1])
for lab, c in zip(grp_prop.index, pal):
    ax.bar(range(4), grp_prop.loc[lab], bottom=bottom, label=lab, color=c,
           width=0.6, edgecolor="white", linewidth=0.3)
    bottom += grp_prop.loc[lab].to_numpy()
ax.set_xticks(range(4)); ax.set_xticklabels(grp_prop.columns, fontsize=9)
ax.set_ylabel("Proportion of microglia"); ax.set_ylim(0, 1)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6.5, frameon=False)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig2c_mg_composition.png", dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig2d 每样本比例统计 ----
tab = pd.crosstab([mg.obs["sample"], mg.obs["group"]], mg.obs["mg_label"])
prop = tab.div(tab.sum(axis=1), axis=0).reset_index()
prop.to_csv(f"{OUT}/mg_state_proportions_per_sample.csv", index=False)
rows = []
for tissue in ["Retina", "OpticNerve"]:
    sub = prop[prop["group"].str.startswith(tissue)]
    for lab in order:
        a = sub.loc[sub["group"].str.endswith("MB"), lab].to_numpy()
        b = sub.loc[sub["group"].str.endswith("Sham"), lab].to_numpy()
        if a.sum() == 0 and b.sum() == 0: continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append({"tissue": tissue, "state": lab, "prop_MB": a.mean(),
                     "prop_Sham": b.mean(),
                     "log2FC": np.log2((a.mean() + 1e-4) / (b.mean() + 1e-4)),
                     "pvalue": p})
st = pd.DataFrame(rows)
st["padj"] = st.groupby("tissue")["pvalue"].transform(
    lambda p: np.minimum(p * len(p) / (np.argsort(np.argsort(p)) + 1), 1))
st.to_csv(f"{OUT}/mg_state_proportion_stats.csv", index=False)
print(st.round(4).to_string(), flush=True)
log("Fig2 完成")

# ---- Fig3a: UMAP by DPT ----
sc.pl.umap(mg, color="dpt_pseudotime", cmap="magma", show=False, frameon=False, size=4)
plt.savefig(f"{OUT}/Fig3a_UMAP_dpt.png", dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig3b: DPT 密度曲线 ----
fig, ax = plt.subplots(figsize=(6.5, 4))
for gname, color in [("Retina-Sham", "#5C988B"), ("Retina-MB", "#C0532F"),
                     ("OpticNerve-Sham", "#8FA8A0"), ("OpticNerve-MB", "#B8860B")]:
    v = mg.obs.loc[mg.obs["group"] == gname, "dpt_pseudotime"].to_numpy()
    v = v[np.isfinite(v)]
    kde = stats.gaussian_kde(v)
    xs = np.linspace(0, np.nanmax(mg.obs["dpt_pseudotime"]), 300)
    ax.plot(xs, kde(xs), color=color, lw=2, label=f"{gname} (n={len(v)})")
ax.set_xlabel("Pseudotime (DPT)"); ax.set_ylabel("Density")
ax.legend(fontsize=8, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig3b_dpt_density.png", dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig3c: 分组箱线图 + 统计 ----
fig, axes = plt.subplots(1, 2, figsize=(9, 4))
stat_rows = []
for ax, tissue in zip(axes, ["Retina", "OpticNerve"]):
    a = mg.obs.loc[(mg.obs["tissue"] == tissue) & (mg.obs["condition"] == "MB"), "dpt_pseudotime"].dropna()
    b = mg.obs.loc[(mg.obs["tissue"] == tissue) & (mg.obs["condition"] == "Sham"), "dpt_pseudotime"].dropna()
    bp = ax.boxplot([b, a], labels=["Sham", "MB"], patch_artist=True, widths=0.5,
                    medianprops=dict(color="#3E2723"))
    for patch, c in zip(bp["boxes"], ["#8FA8A0", "#C0532F"]):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    u, p = stats.mannwhitneyu(a, b)
    d = a.median() - b.median()
    stat_rows.append({"tissue": tissue, "median_MB": a.median(), "median_Sham": b.median(),
                      "delta_median": d, "mwu_p": p})
    ax.set_title(f"{tissue}\nΔmedian={d:.3f}, MWU p={p:.1e}", fontsize=10)
    ax.set_ylabel("Pseudotime")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig3c_dpt_boxplot.png", dpi=300, bbox_inches="tight"); plt.close()
pd.DataFrame(stat_rows).to_csv(f"{OUT}/dpt_group_stats.csv", index=False)
print(pd.DataFrame(stat_rows).round(4).to_string(), flush=True)

# ---- Fig3d: 轨迹基因热图 ----
traj = pd.read_csv(f"{OUT}/traj_gene_correlations.csv")
KNOWN = ["P2ry12", "Tmem119", "Sall1", "Cx3cr1", "Apoe", "Ctsb", "Ctsd", "Cd63",
         "Lyz2", "Trem2", "Spp1", "Cst7", "Cd74", "H2-Aa", "Ifit3", "Isg15",
         "Ccl5", "Lgals3", "Fth1", "Mki67"]
top_neg = traj.head(12)["gene"].tolist()
top_pos = traj.tail(12)["gene"].tolist()[::-1]
genes_hm = list(dict.fromkeys(top_neg + KNOWN + top_pos))
genes_hm = [g for g in genes_hm if g in mg.var_names][:48]

dpt = mg.obs["dpt_pseudotime"].to_numpy()
valid = np.isfinite(dpt)
order_cells = np.argsort(dpt[valid])
Xv = mg[valid][:, genes_hm].X
if sp.issparse(Xv): Xv = Xv.toarray()
Xv = Xv[order_cells]
# 分箱平滑
nb = 120
bins = np.array_split(np.arange(Xv.shape[0]), nb)
Hm = np.vstack([Xv[b].mean(axis=0) for b in bins])
Hm = np.log1p(Hm * 100)  # 视觉压缩
Hm = (Hm - Hm.mean(axis=0)) / (Hm.std(axis=0) + 1e-9)
fig, ax = plt.subplots(figsize=(7, 8))
im = ax.imshow(Hm.T, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
               extent=[0, 1, len(genes_hm), 0])
ax.set_yticks(np.arange(len(genes_hm)) + 0.5)
ax.set_yticklabels(genes_hm, fontsize=6)
ax.set_xlabel("Pseudotime (homeostatic → activated)")
plt.colorbar(im, label="z-scored expr", shrink=0.5)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig3d_traj_heatmap.png", dpi=300, bbox_inches="tight"); plt.close()
log("Fig3 完成")
print(json.dumps({"figs": sorted(os.listdir(OUT))}, ensure_ascii=False), flush=True)
