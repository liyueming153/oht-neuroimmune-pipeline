#!/usr/bin/env python3
"""19_gse304931_module_scores_4groups.py  (added in revision R1)

Per-sample module scores across the four genotype-by-condition groups of
GSE304931 (FACS-sorted microglia; WT, WT+IOP, Spp1-KO, Spp1-KO+IOP; n=4 each),
with baseline (KO vs WT, no IOP) and IOP-effect-by-genotype tests.

Reviewer 2 asked for sample-level module statistics instead of gene-level
summaries only; this script provides them (Supplementary Figure S5).

Inputs
------
data/GSE304931.csv.gz                          normalized count table (Ensembl rows)
results/05_cci/gse304931_ensembl_symbol_map.csv  Ensembl -> symbol map
    (written by 10_spp1ko_bulk.py; if absent, mygene is used to rebuild it)

Outputs
-------
results/05_cci/GSE304931_module_scores_4groups.csv
    long table: module, group, sample, score, n_genes (96 rows)
results/05_cci/GSE304931_baseline_KO_vs_WT.csv
    per-module Welch t-test of KO vs WT at baseline
results/05_cci/GSE304931_IOP_effect_by_genotype.csv
    per-module IOP effect within each genotype (Welch t-test)
results/05_cci/FigS5_GSE304931_module_scores_4groups.png

Module lists are the extended GSE304931 lists documented in Online Resource 1
(Table S1 note); DAM (excl. Spp1) is a sensitivity control for the mechanical
absence of Spp1 transcripts in knockout samples.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

IN_CSV = "data/GSE304931.csv.gz"
MAP_CSV = "results/05_cci/gse304931_ensembl_symbol_map.csv"
OUT = "results/05_cci"

MODULES = {
    "DAM": ["Apoe", "Spp1", "Lgals3", "Cst7", "Lpl", "Ctsb", "Ctsd", "Ctsl",
            "Trem2", "Axl", "Itgax", "Cd9", "Lilrb4", "Clec7a", "Gpnmb",
            "Fabp5", "Lgals1", "Csf1", "B2m", "Fth1"],
    "IFN": ["Ifit3", "Isg15", "Usp18", "Irf7", "Ifit1", "Ifi44", "Oasl2",
            "Rsad2", "Stat1", "Mx1", "Ifit2", "Oas1a", "Bst2", "Ifi27l2a"],
    "MHCII": ["Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1"],
    "Homeostatic": ["P2ry12", "Tmem119", "Sall1", "Hexb", "Selplg", "Cx3cr1",
                    "Siglech", "Mafb", "Csf1r", "Tgfbr1"],
    "Chemokine": ["Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl12", "Cxcl10", "Cxcl2",
                  "Cxcl16"],
}
MODULES["DAM (excl. Spp1)"] = [g for g in MODULES["DAM"] if g != "Spp1"]

GROUPS = {  # display name -> sample columns
    "WT": ["WT_1", "WT_2", "WT_3", "WT_4"],
    "WT+IOP": ["IOP_1", "IOP_2", "IOP_3", "IOP_4"],
    "Spp1-KO": ["Spp1KO_1", "Spp1KO_2", "Spp1KO_3", "Spp1KO_4"],
    "Spp1-KO+IOP": ["Spp1KO_IOP_1", "Spp1KO_IOP_2", "Spp1KO_IOP_3", "Spp1KO_IOP_4"],
}


def load_symbol_map(ensembl_ids):
    if os.path.exists(MAP_CSV):
        m = pd.read_csv(MAP_CSV)
        return dict(zip(m["ensembl"], m["symbol"]))
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(list(ensembl_ids), scopes="ensembl.gene",
                       fields="symbol", species="mouse")
    m = {r["query"]: r.get("symbol") for r in res if "symbol" in r}
    pd.DataFrame({"ensembl": list(m), "symbol": list(m.values())}).to_csv(
        MAP_CSV, index=False)
    return m


def main():
    df = pd.read_csv(IN_CSV)
    df["ensembl"] = df["gene"].str.split(".").str[0]
    sym = load_symbol_map(df["ensembl"].unique())
    df["symbol"] = df["ensembl"].map(sym)
    df = df.dropna(subset=["symbol"])
    cols = [c for gcols in GROUPS.values() for c in gcols]
    expr = np.log2(df.set_index("symbol")[cols].astype(float) + 1)

    rows = []
    for mod, genes in MODULES.items():
        use = [g for g in genes if g in expr.index]
        for gname, gcols in GROUPS.items():
            for c in gcols:
                rows.append({"module": mod, "group": gname, "sample": c,
                             "score": expr.loc[use, c].mean(),
                             "n_genes": len(use)})
    ms = pd.DataFrame(rows)
    ms.to_csv(f"{OUT}/GSE304931_module_scores_4groups.csv", index=False)

    base_rows, iop_rows = [], []
    for mod in MODULES:
        sub = ms[ms.module == mod]
        wt = sub.loc[sub.group == "WT", "score"].to_numpy()
        ko = sub.loc[sub.group == "Spp1-KO", "score"].to_numpy()
        wi = sub.loc[sub.group == "WT+IOP", "score"].to_numpy()
        ki = sub.loc[sub.group == "Spp1-KO+IOP", "score"].to_numpy()
        base_rows.append({"module": mod, "WT": wt.mean(), "KO": ko.mean(),
                          "delta": ko.mean() - wt.mean(),
                          "welch_p": stats.ttest_ind(ko, wt, equal_var=False).pvalue})
        iop_rows.append({"module": mod,
                         "d_WT": wi.mean() - wt.mean(),
                         "p_WT": stats.ttest_ind(wi, wt, equal_var=False).pvalue,
                         "d_Spp1-KO": ki.mean() - ko.mean(),
                         "p_Spp1-KO": stats.ttest_ind(ki, ko, equal_var=False).pvalue})
    pd.DataFrame(base_rows).to_csv(f"{OUT}/GSE304931_baseline_KO_vs_WT.csv",
                                   index=False)
    pd.DataFrame(iop_rows).to_csv(f"{OUT}/GSE304931_IOP_effect_by_genotype.csv",
                                  index=False)

    # ---- Supplementary Figure S5 ----
    order = ["DAM", "DAM (excl. Spp1)", "IFN", "MHCII", "Homeostatic",
             "Chemokine"]
    gorder = list(GROUPS)
    colors = {"WT": "#9e9e9e", "WT+IOP": "#6baed6", "Spp1-KO": "#d9d9d9",
              "Spp1-KO+IOP": "#fd8d3c"}
    fig, axes = plt.subplots(1, 6, figsize=(16, 3.6), sharey=False)
    for ax, mod in zip(axes, order):
        data = [ms.loc[(ms.module == mod) & (ms.group == g), "score"].to_numpy()
                for g in gorder]
        bp = ax.boxplot(data, widths=0.6, patch_artist=True,
                        medianprops=dict(color="black"))
        for patch, g in zip(bp["boxes"], gorder):
            patch.set_facecolor(colors[g])
        for i, vals in enumerate(data):
            ax.scatter(np.full(len(vals), i + 1), vals, c="black", s=22,
                       zorder=3)
        ax.set_title(mod)
        ax.set_xticklabels(["WT", "WT\n+IOP", "KO", "KO\n+IOP"])
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("Module score (mean log2(norm+1))")
    fig.suptitle("GSE304931 sorted microglia: module scores across the four "
                 "genotype-by-condition groups (n=4 per group)")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(f"{OUT}/FigS5_GSE304931_module_scores_4groups.png", dpi=300,
                facecolor="white")
    print("saved module-score tables and Fig S5")


if __name__ == "__main__":
    main()
