#!/usr/bin/env python3
"""20_spp1_receptor_benchmark.py  (added in revision R1)

Benchmark of predicted SPP1 receptors among myeloid-to-lymphocyte edges in the
OHT retina, plus per-method convergence statistics for the leading SPP1-CD44
edges and Spp1/Cd44 expression by cell type and condition.

Reviewer 2 asked why CD44 was prioritized over other annotated SPP1 receptors;
this script documents that comparison (Supplementary Tables S2-S4,
Supplementary Figure S6).

Inputs
------
results/05_cci/liana_all_interactions.csv   written by 09_cci.py
    (103,619 candidate interactions; columns include source, target,
     ligand_complex, receptor_complex, lr_means, cellphone_pvals, expr_prod,
     scaled_weight, lr_logfc, spec_weight, lrscore, specificity_rank,
     magnitude_rank, tissue, condition)
data/GSE293358_final.h5ad                   annotated CD45+ atlas (for Table S4)

Definitions
-----------
Myeloid senders: Homeostatic MG, Activated MG (Ccl5+), DAM-like MG,
    Iron+ MG (Fth1+), Proliferating MG, BAM, MHCII+ MG, MHCII+ Mac/APC, pDC,
    Neutrophil.
Lymphocyte receivers: gdT cell, NK cell, T cell, B cell.
An edge is called OHT-specific when its LIANA consensus specificity_rank
    is <= 0.05 in the microbead (OHT) retina subset.

Outputs
-------
results/05_cci/TableS2_SPP1_receptor_benchmark.csv
results/05_cci/TableS3_SPP1_CD44_per_method.csv
results/05_cci/TableS4_Spp1_Cd44_expression.csv
results/05_cci/FigS6_SPP1_receptor_benchmark.png
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LIANA_CSV = "results/05_cci/liana_all_interactions.csv"
ATLAS_H5AD = "data/GSE293358_final.h5ad"
OUT = "results/05_cci"

RECEPTORS = {  # display name -> LIANA receptor_complex
    "CD44": "Cd44",
    "integrin aVb1": "Itgav_Itgb1",
    "integrin a4b1": "Itga4_Itgb1",
    "integrin aVb5": "Itgav_Itgb5",
    "integrin a9b1": "Itga9_Itgb1",
    "PTGER4": "Ptger4",
    "S1PR1": "S1pr1",
    "integrin aVb3": "Itgav_Itgb3",
    "integrin a5b1": "Itga5_Itgb1",
}
# Greek-letter display names used in the deposited table / figure
RECEPTORS = {k.replace("a", "\u03b1", 1).replace("b", "\u03b2", 1)
             if k.startswith("integrin") else k: v
             for k, v in RECEPTORS.items()}

MYELOID = ["Homeostatic MG", "Activated MG (Ccl5+)", "DAM-like MG",
           "Iron+ MG (Fth1+)", "Proliferating MG", "BAM", "MHCII+ MG",
           "MHCII+ Mac/APC", "pDC", "Neutrophil"]
LYMPH = ["gdT cell", "NK cell", "T cell", "B cell"]
SPEC_CUTOFF = 0.05


def build_table_s2(li):
    sp = li[(li.ligand_complex.str.lower() == "spp1")
            & (li.tissue == "Retina") & (li.condition == "MB")]
    ml = sp[sp.source.isin(MYELOID) & sp.target.isin(LYMPH)]
    rows = []
    for name, rc in RECEPTORS.items():
        e = ml[ml.receptor_complex == rc]
        a = sp[sp.receptor_complex == rc]
        rows.append({
            "SPP1 receptor": name,
            "mean_lr_myeloid_to_lymph": round(e.lr_means.mean(), 3),
            "sem": round(e.lr_means.sem(), 3),
            "n_edges": len(e),
            "n_specific": int((e.specificity_rank <= SPEC_CUTOFF).sum()),
            "all_edges_n": len(a),
            "all_edges_mean_lr": round(a.lr_means.mean(), 3),
        })
    return (pd.DataFrame(rows)
            .sort_values("mean_lr_myeloid_to_lymph", ascending=False)
            .reset_index(drop=True))


def build_table_s3(li, top_n=8):
    sp = li[(li.ligand_complex.str.lower() == "spp1")
            & (li.receptor_complex == "Cd44")
            & (li.tissue == "Retina") & (li.condition == "MB")]
    ml = sp[sp.source.isin(MYELOID) & sp.target.isin(LYMPH)]
    top = ml.sort_values("lr_means", ascending=False).head(top_n)
    return pd.DataFrame({
        "source": top.source, "target": top.target,
        "LR-means (consensus)": top.lr_means.round(3),
        "CellPhoneDB p": top.cellphone_pvals.round(3),
        "NATMI expr_prod": top.expr_prod.round(3),
        "Connectome scaled_weight": top.scaled_weight.round(3),
        "logFC weight": top.lr_logfc.round(3),
        "SingleCellSignalR lrscore": top.lrscore.round(3),
        "specificity_rank": top.specificity_rank.round(3),
        "magnitude_rank": top.magnitude_rank.round(3),
    }).reset_index(drop=True)


def build_table_s4():
    """Retina-only mean log-normalized expression (log1p of CP10K counts) and
    percent expressing, per cell class and condition, for Spp1 and Cd44.
    gdT cells and pDCs have no sham-retina cells and appear for the OHT
    condition only; the MHCII+ Mac/APC sham estimate rests on 5 cells."""
    import scanpy as sc
    ad = sc.read_h5ad(ATLAS_H5AD, backed="r")
    sub = ad[ad.obs["tissue"] == "Retina", ["Spp1", "Cd44"]].to_memory()
    df = sub.obs[["celltype", "condition"]].copy()
    import scipy.sparse as sp_
    X = sub.X.toarray() if sp_.issparse(sub.X) else np.asarray(sub.X)
    df["Spp1"], df["Cd44"] = X[:, 0], X[:, 1]
    rows = []
    for (ct, cond), grp in df.groupby(["celltype", "condition"], observed=True):
        rows.append({"celltype": ct, "condition": cond, "n_cells": len(grp),
                     "Spp1": round(grp.Spp1.mean(), 3),
                     "Spp1_pct": round((grp.Spp1 > 0).mean(), 3),
                     "Cd44": round(grp.Cd44.mean(), 3),
                     "Cd44_pct": round((grp.Cd44 > 0).mean(), 3)})
    ct_order = ["Homeostatic MG", "MHCII+ MG", "Activated MG (Ccl5+)",
                "DAM-like MG", "Iron+ MG (Fth1+)", "Proliferating MG", "BAM",
                "MHCII+ Mac/APC", "pDC", "T cell", "gdT cell", "NK cell",
                "B cell", "Neutrophil"]
    out = pd.DataFrame(rows)
    out["celltype"] = pd.Categorical(out["celltype"], categories=ct_order,
                                     ordered=True)
    out["condition"] = pd.Categorical(out["condition"],
                                      categories=["Sham", "MB"], ordered=True)
    return out.sort_values(["celltype", "condition"]).reset_index(drop=True)


def fig_s6(ts2):
    d = ts2.iloc[::-1]  # horizontal bars, CD44 on top after sorting desc
    colors = ["#d95f02" if r == "CD44" else "#8dd3c7" for r in d["SPP1 receptor"]]
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    ax.barh(d["SPP1 receptor"], d.mean_lr_myeloid_to_lymph,
            xerr=d.sem, color=colors, edgecolor="white",
            error_kw=dict(ecolor="0.35", capsize=3, lw=1))
    for i, (_, r) in enumerate(d.iterrows()):
        ax.text(0.70, i, f"{r.n_specific}/{r.n_edges}", va="center",
                fontsize=9)  # annotation column, clear of the error bars
    ax.text(0.70, len(d) - 0.2, "specific/\nedges", fontsize=8, ha="left",
            va="bottom", color="0.3")
    ax.set_xlim(0, 0.78)
    ax.set_xlabel("Mean LR-means (myeloid -> lymphocyte edges, retina OHT)")
    ax.set_title("Benchmarking of predicted SPP1 receptors among "
                 "myeloid->lymphocyte edges\n(GSE293358 retina, OHT, CD44 in red)")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/FigS6_SPP1_receptor_benchmark.png", dpi=300,
                facecolor="white")


def main():
    li = pd.read_csv(LIANA_CSV)
    ts2 = build_table_s2(li)
    ts2.to_csv(f"{OUT}/TableS2_SPP1_receptor_benchmark.csv", index=False)
    build_table_s3(li).to_csv(f"{OUT}/TableS3_SPP1_CD44_per_method.csv",
                              index=False)
    build_table_s4().to_csv(f"{OUT}/TableS4_Spp1_Cd44_expression.csv",
                            index=False)
    fig_s6(ts2)
    print("saved Tables S2-S4 and Fig S6")


if __name__ == "__main__":
    main()
