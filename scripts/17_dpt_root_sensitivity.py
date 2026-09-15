#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""阶段3-补充：DPT 根细胞选择敏感性分析（15 个替代根）

设计：在稳态小胶质亚群内挑选 15 个替代根细胞 ——
  (1) 稳态得分最高的 10 个细胞（top10_homeostatic_score）；
  (2) 稳态细胞内随机抽取 5 个（random_homeostatic, default_rng(seed=0)）。
稳态亚群 = mg_state_signatures.csv 中 Homeostatic 得分 top4 的状态（2/1/0/6）。
稳态得分 = 10 个稳态标志基因（P2ry12, Tmem119, Cx3cr1, Sall1, Siglech,
Hexb, Fcrls, Olfml3, Gpr34, P2ry13）在 log-normalized X 上的均值。
逐根重置 uns["iroot"] 并重跑 sc.tl.dpt（复用已有 neighbors / X_diffmap），
与参考 dpt_pseudotime 计算 Spearman 秩相关。

输入：data/GSE293358_mg.h5ad
      （含 X_pca / X_diffmap / neighbors / iroot / dpt_pseudotime / mg_state；
        由 05_mg_trajectory.py 生成；分卷合并：cat GSE293358_mg.h5ad.part_* > GSE293358_mg.h5ad）
      results/03_mg/mg_state_signatures.csv
输出：results/03_mg/dpt_root_sensitivity.csv
      （手稿对应运行结果：median rho = 0.44, range 0.22-0.57）

用法：python3 scripts/17_dpt_root_sensitivity.py [mg.h5ad] [signatures.csv] [out.csv]
"""
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

H5AD = sys.argv[1] if len(sys.argv) > 1 else "data/GSE293358_mg.h5ad"
SIG  = sys.argv[2] if len(sys.argv) > 2 else "results/03_mg/mg_state_signatures.csv"
OUT  = sys.argv[3] if len(sys.argv) > 3 else "results/03_mg/dpt_root_sensitivity.csv"

HOMEO_GENES = ["P2ry12", "Tmem119", "Cx3cr1", "Sall1", "Siglech",
               "Hexb", "Fcrls", "Olfml3", "Gpr34", "P2ry13"]

m = sc.read_h5ad(H5AD)
reference_dpt = m.obs["dpt_pseudotime"].to_numpy().copy()

# 稳态亚群：签名表 Homeostatic 得分 top4 的 mg_state
sig = pd.read_csv(SIG)
top4 = (sig.assign(mg_state=sig["mg_state"].astype(str))
           .sort_values("Homeostatic", ascending=False)["mg_state"]
           .head(4).tolist())
homeo_idx = np.where(m.obs["mg_state"].astype(str).isin(top4).to_numpy())[0]

# 稳态得分（log-normalized X 上 10 基因均值）
gidx = [m.var_names.get_loc(g) for g in HOMEO_GENES]
hscore = np.asarray(m.X[:, gidx].mean(axis=1)).ravel()

# 15 个替代根：top10 稳态得分 + 5 个随机稳态细胞（seed=0）
top10 = homeo_idx[np.argsort(-hscore[homeo_idx])[:10]]
rand5 = np.random.default_rng(0).choice(homeo_idx, 5, replace=False)
roots = np.concatenate([top10, rand5])
selection = ["top10_homeostatic_score"] * 10 + ["random_homeostatic"] * 5

rows = []
for i, ci in enumerate(roots):
    m.uns["iroot"] = int(ci)
    sc.tl.dpt(m)  # 复用已有 neighbors / X_diffmap
    rho = spearmanr(m.obs["dpt_pseudotime"].to_numpy(), reference_dpt).statistic
    rows.append({"root_id": i + 1,
                 "cell_index": int(ci),
                 "selection": selection[i],
                 "homeostatic_score": float(hscore[ci]),
                 "spearman_rho_vs_reference": float(rho)})
    print(f"root {i + 1:2d}/15  cell {ci:6d}  rho = {rho:.3f}")

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(f"\nmedian rho = {out['spearman_rho_vs_reference'].median():.2f}, "
      f"range {out['spearman_rho_vs_reference'].min():.2f}-{out['spearman_rho_vs_reference'].max():.2f}")
print("saved ->", OUT)
