#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""阶段4-Part1：GSE241782 时间序列验证（列映射校验 + 胶质反应性评分）"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/tmp/results/04_glia"
os.makedirs(OUT, exist_ok=True)

# GSM 顺序 -> 列名（组织内按 GSM 升序）
META = {}
uon = ["2W_Glaucoma_M","2W_Crush_M","2W_Crush_F","2W_Glaucoma_F","3D_Crush_F","3D_Crush_M",
       "3D_Glaucoma_M","3D_Glaucoma_F","6W_Glaucoma_F","6W_Glaucoma_M","0D_Naive_F","0D_Naive_M"]
mon = ["2W_Crush_M","2W_Crush_F","2W_Glaucoma_M","2W_Glaucoma_F","3D_Crush_F","3D_Crush_M",
       "3D_Glaucoma_M","3D_Glaucoma_F","6W_Glaucoma_M","6W_Glaucoma_F","0D_Naive_F","0D_Naive_M"]
ret = ["2W_Crush_F","2W_Crush_M","2W_Glaucoma_M","2W_Glaucoma_F","3D_Crush_M","3D_Crush_F",
       "3D_Glaucoma_M","3D_Glaucoma_F","6W_Glaucoma_M","6W_Glaucoma_F","0D_Naive_F","0D_Naive_M"]
for i in range(12):
    META[f"UON-{i+1:02d}"] = uon[i]; META[f"MON-{i+1:02d}"] = mon[i]; META[f"R-{i+1:02d}"] = ret[i]

df = pd.read_csv("/mnt/agents/output/oht_project/data/GSE241782_gene_count_matrix.csv.gz", index_col=0)
print("矩阵:", df.shape, flush=True)
# 基因名索引处理：ENSMUSG|Symbol -> Symbol
df.index = df.index.astype(str).str.split("|").str[-1]
df = df.groupby(df.index).max()

# log-CPM
lib = df.sum(axis=0)
cpm = df.div(lib, axis=1) * 1e6
logcpm = np.log2(cpm + 1)

# ---- Xist 校验性别 ----
xist = logcpm.loc["Xist"] if "Xist" in logcpm.index else None
if xist is not None:
    sex_df = pd.DataFrame({"col": xist.index, "Xist": xist.values})
    sex_df["sex_pred"] = np.where(sex_df["Xist"] > 2, "F", "M")
    sex_df["sex_meta"] = [META[c].split("_")[-1] for c in sex_df["col"]]
    sex_df["match"] = sex_df["sex_pred"] == sex_df["sex_meta"]
    print(sex_df.to_string(index=False), flush=True)
    print(f"性别校验匹配率: {sex_df['match'].mean()*100:.0f}%", flush=True)

meta_df = pd.DataFrame({"col": list(META.keys()),
    "tissue": [c.split("-")[0].replace("R", "Retina") for c in META],
    "time": [META[c].split("_")[0] for c in META],
    "injury": [META[c].split("_")[1] for c in META],
    "sex": [META[c].split("_")[2] for c in META]})
meta_df.to_csv(f"{OUT}/GSE241782_sample_meta.csv", index=False)

# ---- 评分签名 ----
SIG = {
    "Pan_reactive": ["Gfap", "Vim", "Serpina3n", "Clu", "S100a6", "Apoe"],
    "A1_neurotoxic": ["H2-D1", "Serping1", "Ggta1", "Iigp1", "Gbp2", "Fbln5",
                      "Fkbp5", "Psmb8", "Srgn", "Amigo2"],
    "A2_neuroprotective": ["Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd14",
                           "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd109", "Ptgs2"],
    "Complement": ["C1qa", "C1qb", "C1qc", "C3", "C4b"],
    "DAM_MG": ["Apoe", "Ctsb", "Ctsd", "Cd63", "Lyz2", "Lgals3", "Spp1", "Cst7"],
    "IFN_module": ["Ifit3", "Isg15", "Irf7", "Stat1", "Usp18", "Ifit1"],
    "T_NK_infil": ["Cd3e", "Cd8a", "Nkg7", "Gzma"],
    "Neutrophil": ["S100a8", "S100a9", "Mpo"],
}
scores = pd.DataFrame(index=logcpm.columns)
for sname, genes in SIG.items():
    g_use = [g for g in genes if g in logcpm.index]
    scores[sname] = logcpm.loc[g_use].mean(axis=0)
    print(sname, "genes found:", len(g_use), "/", len(genes), flush=True)
scores.index.name = "col"
sc = scores.reset_index().merge(meta_df, on="col")
sc.to_csv(f"{OUT}/GSE241782_signature_scores.csv", index=False)

# ---- Fig4a: 时间序列（青光眼组） ----
fig, axes = plt.subplots(2, 4, figsize=(16, 7))
plot_sigs = ["Pan_reactive", "A1_neurotoxic", "A2_neuroprotective", "Complement",
             "DAM_MG", "IFN_module", "T_NK_infil", "Neutrophil"]
tissues = [("Retina", "#205B57"), ("UON", "#B8860B"), ("MON", "#6B3A5B")]
times = ["0D", "3D", "2W", "6W"]
for ax, sname in zip(axes.ravel(), plot_sigs):
    for tissue, color in tissues:
        sub = sc[(sc["tissue"] == tissue) & (sc["injury"].isin(["Naive", "Glaucoma"]))]
        sub = sub.copy()
        sub["time"] = pd.Categorical(sub["time"], categories=times, ordered=True)
        gm = sub.groupby("time", observed=True)[sname].agg(["mean", "sem"])
        ax.errorbar(range(4), gm["mean"], yerr=gm["sem"], marker="o", ms=4,
                    color=color, label=tissue, capsize=3, lw=1.5)
    ax.set_xticks(range(4)); ax.set_xticklabels(times)
    ax.set_title(sname, fontsize=11)
    ax.set_xlabel("Time after glaucoma induction")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes.ravel()[0].legend(fontsize=8, frameon=False)
plt.tight_layout()
plt.savefig(f"{OUT}/Fig4a_timecourse_scores.png", dpi=300, bbox_inches="tight"); plt.close()

# 青光眼 3D vs naive 的倍数统计
rows = []
for tissue in ["Retina", "UON", "MON"]:
    for sname in plot_sigs:
        a = sc[(sc["tissue"] == tissue) & (sc["time"] == "3D") & (sc["injury"] == "Glaucoma")][sname]
        b = sc[(sc["tissue"] == tissue) & (sc["injury"] == "Naive")][sname]
        if len(a) and len(b):
            t, p = stats.ttest_ind(a, b, equal_var=False)
            rows.append({"tissue": tissue, "signature": sname,
                         "mean_3D": a.mean(), "mean_naive": b.mean(),
                         "delta": a.mean() - b.mean(), "pvalue": p})
st = pd.DataFrame(rows)
st.to_csv(f"{OUT}/GSE241782_3D_vs_naive_stats.csv", index=False)
print(st.round(3).to_string(index=False), flush=True)
print("PART1 完成", flush=True)
