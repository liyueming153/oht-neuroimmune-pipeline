#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""阶段4-Part2：GSE306785 分选 Müller 细胞 OHT vs SHAM 差异与反应性"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/tmp/results/04_glia"
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv("/mnt/agents/output/oht_project/data/GSE306785_normalized_macs.csv.gz", index_col=0)
df.index = df.index.astype(str)
df = df.groupby(df.index).max()
print("Müller 矩阵:", df.shape, df.columns.tolist(), flush=True)
print(df.describe().loc[["mean", "max"]].round(1).to_string(), flush=True)

SHAM = ["SHAM_1", "SHAM_2", "SHAM_3"]
OHT = ["OHT_1", "OHT_2", "OHT_3"]
logdf = np.log2(df + 1)

# ---- 差异分析（log 空间 Welch t + log2FC on normalized counts） ----
res = []
for g in df.index:
    a = logdf.loc[g, OHT].to_numpy()
    b = logdf.loc[g, SHAM].to_numpy()
    if df.loc[g].max() < 10:   # 过滤低表达
        continue
    t, p = stats.ttest_ind(a, b, equal_var=False)
    lfc = np.log2((df.loc[g, OHT].mean() + 1) / (df.loc[g, SHAM].mean() + 1))
    res.append({"gene": g, "log2FC": lfc, "pvalue": p,
                "mean_OHT": df.loc[g, OHT].mean(), "mean_SHAM": df.loc[g, SHAM].mean()})
de = pd.DataFrame(res).sort_values("pvalue")
de["padj"] = np.minimum(de["pvalue"] * len(de) /
                        (pd.Series(de["pvalue"].rank()).to_numpy()), 1)
de.to_csv(f"{OUT}/GSE306785_Muller_DE.csv", index=False)
nsig = ((de["padj"] < 0.05) & (de["log2FC"].abs() > 0.5)).sum()
print(f"显著差异基因 (FDR<0.05 & |log2FC|>0.5): {nsig}", flush=True)
print(de.head(20).round(3).to_string(index=False), flush=True)

# ---- 反应性评分 ----
SIG = {
    "Pan_reactive": ["Gfap", "Vim", "Serpina3n", "Clu", "S100a6", "Apoe"],
    "A1_neurotoxic": ["H2-D1", "Serping1", "Ggta1", "Iigp1", "Gbp2", "Fbln5",
                      "Fkbp5", "Psmb8", "Srgn", "Amigo2"],
    "A2_neuroprotective": ["Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd14",
                           "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd109", "Ptgs2"],
    "Complement": ["C1qa", "C1qb", "C1qc", "C3", "C4b"],
    "Chemokine_recruit": ["Ccl2", "Ccl5", "Cxcl10", "Cxcl12", "Cxcl16", "Il6", "Il1b", "Tnf"],
}
sc = {}
for sname, genes in SIG.items():
    g_use = [g for g in genes if g in df.index]
    sc[sname] = logdf.loc[g_use].mean(axis=0)
sc = pd.DataFrame(sc)
sc["cond"] = ["OHT"] * 3 + ["SHAM"] * 3
sc.to_csv(f"{OUT}/GSE306785_Muller_scores.csv")

fig, ax = plt.subplots(figsize=(8, 4.5))
x = np.arange(len(SIG))
for i, (cond, color, off) in enumerate([("SHAM", "#8FA8A0", -0.18), ("OHT", "#C0532F", 0.18)]):
    m = [sc.loc[sc["cond"] == cond, s].mean() for s in SIG]
    e = [sc.loc[sc["cond"] == cond, s].sem() for s in SIG]
    ax.bar(x + off, m, 0.35, yerr=e, color=color, label=cond, capsize=3, alpha=0.9)
# 显著性标注
for j, sname in enumerate(SIG):
    a = sc.loc[sc["cond"] == "OHT", sname]; b = sc.loc[sc["cond"] == "SHAM", sname]
    t, p = stats.ttest_ind(a, b, equal_var=False)
    if p < 0.05:
        ax.text(j, max(a.mean(), b.mean()) + 0.15, f"p={p:.3f}", ha="center", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(list(SIG.keys()), rotation=20, ha="right", fontsize=8)
ax.set_ylabel("Mean log2(norm+1)")
ax.legend(frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig4b_Muller_scores.png", dpi=300, bbox_inches="tight"); plt.close()

# ---- Fig4c: 火山图 ----
fig, ax = plt.subplots(figsize=(6, 5.5))
sig = (de["padj"] < 0.05) & (de["log2FC"].abs() > 0.5)
ax.scatter(de.loc[~sig, "log2FC"], -np.log10(de.loc[~sig, "pvalue"]),
           s=3, c="#BBBBBB", rasterized=True)
ax.scatter(de.loc[sig, "log2FC"], -np.log10(de.loc[sig, "pvalue"]),
           s=5, c="#C0532F", rasterized=True)
KEY = ["Gfap", "Vim", "Apoe", "Clu", "C3", "Ccl2", "Cxcl10", "Spp1", "Serpina3n",
       "S100a6", "Lcn2", "Ednrb", "Aqp4", "Rlbp1", "Glul", "Il6st", "Osmr", "Stat3",
       "Cx3cl1", "Tgfb1", "Lif", "Socs3"]
for g in KEY:
    if g in de["gene"].values:
        r = de[de["gene"] == g].iloc[0]
        ax.annotate(g, (r["log2FC"], -np.log10(r["pvalue"])), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
ax.axvline(0.5, ls="--", lw=0.7, c="gray"); ax.axvline(-0.5, ls="--", lw=0.7, c="gray")
ax.axhline(-np.log10(0.05), ls="--", lw=0.7, c="gray")
ax.set_xlabel("log2FC (OHT vs SHAM, sorted Müller cells)")
ax.set_ylabel("-log10 p")
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig4c_Muller_volcano.png", dpi=300, bbox_inches="tight"); plt.close()

# 关键基因表
key_de = de[de["gene"].isin(KEY)].copy()
key_de.to_csv(f"{OUT}/GSE306785_keygenes.csv", index=False)
print(key_de.round(3).to_string(index=False), flush=True)
print("PART2 完成", flush=True)
