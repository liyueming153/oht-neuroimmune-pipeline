# OHT Neuroimmune Pipeline

Analysis code for the manuscript:

**"A Multi-Cohort Transcriptomic Re-Analysis of Microglial States and Candidate Immune Signaling in Ocular Hypertension"**

All analyses use public data. No raw data are redistributed here; download instructions are given below.

## Pipeline (run in order)

| Script | Step | Input | Key output |
|---|---|---|---|
| `01_load_qc.py` | Load 13 GSE293358 10x samples, per-sample QC (genes ≥200, mt <10%, sample-specific 5×MAD upper bound), Scrublet doublet calls (expected rate 0.006, threshold 0.25; CD45+-sorted data), sparse merge | `data/GSE293358/` | `out/merged_qc.h5ad`, `out/qc_table.csv` |
| `02_annotate.py` | Normalize (1e4 + log1p), 2,500 HVGs (seurat, batch=sample), PCA(50), k=15 graph, UMAP, Leiden 1.0; marker z-score argmax annotation + raw-expression lineage gates + cell-level lymphoid refinement; contaminant gating | `out/merged_qc.h5ad` | `out/annotated.h5ad` |
| `03_liana.py` | LIANA 1.8.1 `rank_aggregate` (mouseconsensus resource; expr_prop=0.1, min_cells=10, **1,000 permutations**, seed=1337), run per tissue × condition, one subset per process | `out/annotated.h5ad` | `out/liana_<tissue>_<cond>_1000perms.csv.gz` |
| `04_fig5a.py` | Figure 5a: Spp1→CD44 sender×receiver heatmaps (Sham vs MB, retina); cells with specificity_rank > 0.05 masked | `out/liana_Retina_*_1000perms.csv.gz` | `out/fig5a_1000perms.png` |
| `05_bulk_gse241782.py` | Bulk time-course module replication (GSE241782) | `data/GSE241782/` | `out/gse241782_module_scores.csv`, `out/gse241782_module_tests.csv` |
| `06_bulk_gse306785.py` | Sorted Müller glia / microglia CX3CL1–CX3CR1 module analysis (GSE306785) | `data/GSE306785/` | module tests |
| `07_bulk_gse304931.py` | αRGC regenerative cohort checks (GSE304931) | `data/GSE304931/` | module tests |
| `08_pseudotime.py` | Microglia diffusion pseudotime (homeostatic→DAM→IFN/MHCII continuum), gene–pseudotime correlations | `out/annotated.h5ad` | `out/dpt_trajectory_genes.csv`, `out/microglia_dpt.h5ad` |
| `09_mr.py` | Mendelian randomization (GTEx v8 whole-blood eQTL → FinnGen R10 glaucoma) | `data/MR/` | `out/mr_results.csv` |
| `10_lincs.py` | LINCS L1000 signature-reversal mining via Enrichr | internet | `out/lincs_reversal_top30.csv` |
| `fig_sc_all.py` | Figures 1–3 + 5b (publication style: Arial 7 pt, Okabe–Ito palette) + composition/state statistics | `out/annotated.h5ad`, `out/microglia_dpt.h5ad` | `out/fig1a–3d_new.png`, `out/fig5b_new.png`, `out/composition_stats.csv`, `out/mgstate_stats.csv`, `out/dpt_group_tests.csv` |
| `fig_bulk.py`, `fig_bulk2.py` | Figure 4 (bulk module scores, Müller volcano) | bulk result CSVs | `out/fig4a–4c_new.png` |
| `fig5a2.py`, `fig5cd.py`, `fig67.py` | Figures 5a (sender×receiver heatmaps), 5c–d, 6 (MR forest), 7 (LINCS) | results CSVs | `results/figures/fig5a,5c,5d,6,7_new.png` |
| `root_sens.py` | Pseudotime root-choice sensitivity (15 alternative roots, Spearman ρ) | `out/microglia_dpt.h5ad` | `out/root_sensitivity.csv` |

## Reproduce

```bash
pip install -r requirements.txt   # scanpy==1.12.3, liana==1.8.1, scrublet, ...

# 1. download GSE293358 (13 samples x barcodes/features/matrix .gz) into data/GSE293358/
#    from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293358
python3 scripts/01_load_qc.py
python3 scripts/02_annotate.py

# 2. cell–cell communication (one subset per process keeps RAM < 4 GB)
for tc in "Retina Sham" "Retina MB" "OpticNerve Sham" "OpticNerve MB"; do
  set -- $tc
  python3 scripts/03_liana.py "$1" "$2" 1000
done
python3 scripts/04_fig5a.py
```

Paths are configurable via environment variables: `GSE293358_DIR`, `GSE241782_DIR`,
`GSE306785_DIR`, `GSE304931_DIR`, `MR_DIR`, `OUT_DIR` (defaults: `data/<Cohort>`, `out`).

## Key 1,000-permutation LIANA result

Spp1–CD44 is **entirely absent in sham retina** (0 of 12,973 candidates) and becomes
pervasive under microbead-induced OHT (33 edges, 18 passing specificity_rank ≤ 0.05):
neutrophils and activated myeloid cells (DAM-like microglia, MHCII+ macrophages/APCs)
are the senders; γδ T, NK, T, B cells and monocytes are the receivers.
Full per-subset tables are in `results/liana_*_1000perms.csv.gz`.

## Data availability

- GSE293358, GSE241782, GSE306785, GSE304931 — NCBI GEO
- GTEx v8 — gtexportal.org; FinnGen R10 — finngen.fi
- LINCS L1000 — via Enrichr (maayanlab.cloud)

## License

MIT (code). Public datasets retain their original terms of use.
