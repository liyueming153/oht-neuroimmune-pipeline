#!/usr/bin/env python3
"""18_tissue_stratified_dpt.py  (added in revision R1)

Tissue-stratified diffusion-pseudotime (DPT) sensitivity analysis.

Reviewer 2 asked whether the microglial DPT continuum could be confounded by
tissue of origin, because the discovery run pooled retina and optic-nerve
microglia. This script recomputes DPT independently within each tissue with
the identical pipeline (HVG selection, neighbor graph, homeostatic rooting)
and tests whether the state ordering is preserved.

Inputs
------
data/GSE293358_mg.h5ad   microglia subset written by 05_mg_trajectory.py
                         (obs must contain: tissue ("Retina"/"OpticNerve"),
                          sample, condition, mg_state)

Outputs
-------
results/03_mg/dpt_tissue_stratified_state_medians.csv
    Median tissue-specific DPT of each of the 12 microglial states, per tissue.
results/03_mg/FigS4_tissue_stratified_DPT.png
    Supplementary Figure S4 (six panels).

Pipeline (identical in both tissues): normalize_total(1e4) -> log1p ->
module scores (Table S1 lists, computed before HVG selection) ->
2,000 HVGs (sc.pp.highly_variable_genes, flavor="seurat", batch_key="sample")
-> subset to HVGs -> scale (float32) -> PCA (40 components) ->
neighbors (n_neighbors=15, n_pcs=25) -> diffmap (12 components) ->
iroot = cell with the highest homeostatic module score -> DPT ->
min-max scaling to [0, 1] -> UMAP.

Note: the OHT-only subset was also evaluated; in retina it degenerated
(too few cells for a stable graph), so only the all-cell tissue-stratified
trajectories are reported.
"""

import gc

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

IN_H5AD = "data/GSE293358_mg.h5ad"
OUT = "results/03_mg"

# Module gene lists (compact single-cell definitions; Online Resource 1, Table S1)
MODULES = {
    "Homeo": ["P2ry12", "Tmem119", "Cx3cr1", "Sall1", "Siglech",
              "Hexb", "Fcrls", "Olfml3", "Gpr34", "P2ry13"],
    "DAM":   ["Apoe", "Spp1", "Lgals3", "Cst7", "Trem2", "Tyrobp",
              "Gpnmb", "Itgax", "Cd9", "Ctsb", "Ctsd", "Csf1"],
    "IFN":   ["Ifit3", "Isg15", "Usp18", "Ifit1", "Oasl2", "Irf7", "Stat1"],
    "MHCII": ["Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Ciita"],
}


def run_dpt(adata):
    """Memory-safe DPT on one tissue subset; returns per-cell obs DataFrame."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # module scores BEFORE HVG subsetting (mean log-normalized expression)
    for name, genes in MODULES.items():
        use = [g for g in genes if g in adata.var_names]
        adata.obs[f"mod_{name}"] = np.asarray(
            adata[:, use].X.mean(axis=1)).ravel()
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat",
                                batch_key="sample")
    adata = adata[:, adata.var.highly_variable].copy()   # dense only after HVG
    adata.X = adata.X.astype(np.float32)
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=40)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=25)
    sc.tl.diffmap(adata, n_comps=12)
    adata.uns["iroot"] = int(np.argmax(adata.obs["mod_Homeo"].to_numpy()))
    sc.tl.dpt(adata)
    dpt = adata.obs["dpt_pseudotime"].to_numpy()
    adata.obs["dpt"] = (dpt - np.nanmin(dpt)) / (np.nanmax(dpt) - np.nanmin(dpt))
    sc.tl.umap(adata)
    cols = ["tissue", "condition", "mg_state", "dpt"] + [f"mod_{m}" for m in MODULES]
    obs = adata.obs[cols].copy()
    obs["UMAP1"] = adata.obsm["X_umap"][:, 0]
    obs["UMAP2"] = adata.obsm["X_umap"][:, 1]
    return obs


def main():
    adata = sc.read_h5ad(IN_H5AD, backed="r")
    per_tissue, medians, mod_rhos = {}, {}, {}
    for tissue in ["Retina", "OpticNerve"]:
        sub = adata[adata.obs["tissue"] == tissue].to_memory()
        obs = run_dpt(sub)
        per_tissue[tissue] = obs
        medians[tissue] = obs.groupby("mg_state")["dpt"].median()
        mod_rhos[tissue] = {m: spearmanr(obs["dpt"], obs[f"mod_{m}"])[0]
                            for m in MODULES}
        del sub
        gc.collect()
    med = pd.DataFrame({"mg_state": medians["Retina"].index,
                        "Retina": medians["Retina"].values,
                        "OpticNerve": medians["OpticNerve"].values})
    med.to_csv(f"{OUT}/dpt_tissue_stratified_state_medians.csv", index=False)
    rho, p = spearmanr(med["Retina"], med["OpticNerve"])
    print(f"state-median concordance: Spearman rho={rho:.3f}, p={p:.3g}")

    # ---- Supplementary Figure S4 ----
    fig = plt.figure(figsize=(16, 9))
    n = {t: len(per_tissue[t]) for t in per_tissue}
    fig.suptitle("Supplementary Figure S4. Tissue-stratified diffusion-pseudotime "
                 f"sensitivity analysis (retina n={n['Retina']:,} microglia; "
                 f"optic nerve n={n['OpticNerve']:,})")
    panels = [("A", "Retina", "state", "Retina: 12 microglial states"),
              ("C", "Retina", "dpt", "Retina: tissue-specific diffusion pseudotime"),
              ("B", "OpticNerve", "state", "Optic nerve: 12 microglial states"),
              ("D", "OpticNerve", "dpt", "Optic nerve: tissue-specific diffusion pseudotime")]
    pos = {"A": (0, 0), "C": (0, 1), "B": (1, 0), "D": (1, 1)}
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.15, 1.0])
    for letter, tissue, color, title in panels:
        ax = fig.add_subplot(gs[pos[letter]])
        o = per_tissue[tissue]
        if color == "state":
            sc_c = ax.scatter(o.UMAP1, o.UMAP2, c=o.mg_state.astype(int),
                              cmap="tab20", s=2)
            fig.colorbar(sc_c, ax=ax, label="state", ticks=range(12))
        else:
            sc_c = ax.scatter(o.UMAP1, o.UMAP2, c=o.dpt, cmap="viridis",
                              s=2, vmin=0, vmax=1)
            fig.colorbar(sc_c, ax=ax, label="DPT")
        ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
    ax = fig.add_subplot(gs[0, 2])
    ax.scatter(med.Retina, med.OpticNerve, c="tab:blue")
    for _, r in med.iterrows():
        ax.annotate(int(r.mg_state), (r.Retina, r.OpticNerve),
                    xytext=(4, 4), textcoords="offset points")
    lim = [0, 0.25]
    ax.plot([0, 0.25], [0, 0.25], ls="--", c="gray")   # y = x within axis limits
    ax.set_xlim(lim); ax.set_ylim([0, 0.9])
    ax.set_xlabel("state-median DPT, retina-specific trajectory")
    ax.set_ylabel("state-median DPT, optic-nerve-specific trajectory")
    ax.set_title(f"E  state-median concordance\n(Spearman rho={rho:.2f}, "
                 f"p={p:.1e}; dashed: y=x)", loc="left", fontweight="bold")
    ax = fig.add_subplot(gs[1, 2])
    xs = [mod_rhos["Retina"][m] for m in MODULES]
    ys = [mod_rhos["OpticNerve"][m] for m in MODULES]
    ax.scatter(xs, ys, c="tab:orange")
    for m, x, y in zip(MODULES, xs, ys):
        ax.annotate(m, (x, y), xytext=(6, 6), textcoords="offset points")
    ax.plot([-0.4, 0.4], [-0.4, 0.4], ls="--", c="gray")
    ax.axhline(0, c="gray", lw=0.5); ax.axvline(0, c="gray", lw=0.5)
    ax.set_xlim([-0.4, 0.4]); ax.set_ylim([-0.4, 0.4])
    ax.set_xlabel("retina: Spearman rho (module score vs DPT)")
    ax.set_ylabel("optic nerve: Spearman rho (module score vs DPT)")
    ax.set_title("F  module-pseudotime correlation,\nboth tissues (dashed: y=x)",
                 loc="left", fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{OUT}/FigS4_tissue_stratified_DPT.png", dpi=300,
                facecolor="white")
    print("saved", f"{OUT}/FigS4_tissue_stratified_DPT.png")


if __name__ == "__main__":
    main()
