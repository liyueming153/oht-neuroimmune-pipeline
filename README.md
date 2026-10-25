# oht-neuroimmune-pipeline

Analysis code for the manuscript:

**"A Hypothesis-Generating Integration of Single-Cell and Bulk
Transcriptomics Nominates SPP1–CD44 Innate–Adaptive Crosstalk in Ocular
Hypertension"** (revision R1, Functional & Integrative Genomics)

All analyses use public data (NCBI GEO, GTEx, FinnGen, LINCS L1000). No raw data
are redistributed here; download sources are listed under **Data availability**.

## Repository layout

```
scripts/    analysis pipeline, run in numeric order (see below)
results/    per-stage result tables (CSV) and figure panels (PNG)
results/plot_data/   source data behind every figure panel (CSV)
```

## Pipeline (run in order)

| Script | Step | Key output |
| --- | --- | --- |
| `01_qc.py` | Load 13 GSE293358 10x samples, per-sample QC metrics & filtering, Scrublet doublet calls, sparse merge | `results/01_qc/`, QC h5ad |
| `02_cluster.py` | Normalize (1e4 + log1p), HVGs, PCA, kNN graph, UMAP, Leiden, marker scoring | clustered h5ad |
| `03_immune_recluster.py` | Remove contaminant clusters, immune reclustering, signature scoring | immune h5ad |
| `04_annotate_fig1.py` | Final annotation of 14 CD45+ immune classes; Fig. 1 atlas panels | `results/02_cluster/` |
| `05_mg_trajectory.py` | Microglia reclustering into 12 states; diffusion pseudotime (DPT); trajectory gene dynamics | `results/03_mg/`, MG h5ad |
| `06_mg_figures.py` | Fig. 2 (microglial states) and Fig. 3 (pseudotime) panels | `results/03_mg/` |
| `07_glia_bulk.py` | GSE241782 bulk time-course module validation (3D–6W; retina/UON/MON) | `results/04_glia/` |
| `08_muller_bulk.py` | GSE306785 FACS-sorted Müller glia, OHT vs sham | `results/04_glia/` |
| `09_cci.py` | LIANA cell–cell communication (mouseconsensus resource), per tissue × condition | `results/05_cci/` |
| `09b_extract_axes.py` | Focused LIANA axes + key-gene expression by cell type | `results/05_cci/` |
| `10_spp1ko_bulk.py` | GSE304931 Spp1-knockout microglia bulk contrasts | `results/05_cci/` |
| `11_fig5.py` | Fig. 5 panels (from CSVs only) | `results/05_cci/` |
| `12_drugrep.py` | LINCS L1000 signature-reversal mining via Enrichr | `results/06_mr_drug/` |
| `13_mr_instruments.py` | GTEx cis-eQTL instrument extraction for core genes | `results/06_mr_drug/` |
| `14_mr_analysis.py` / `14b_mr_from_hits.py` | Two-sample MR: GTEx retina eQTL (exposure) → FinnGen R10 glaucoma (outcome) | `results/06_mr_drug/` |
| `15_fig67.py` | Fig. 6 (MR forest) and Fig. 7 (drug reversal) panels | `results/06_mr_drug/` |
| `16_export_plot_data.py` | Export per-panel source data (GraphPad-ready) | `results/plot_data/` |
| `17_dpt_root_sensitivity.py` | DPT root-choice sensitivity: 15 alternative roots, Spearman ρ vs reference | `results/03_mg/dpt_root_sensitivity.csv` |
| `18_tissue_stratified_dpt.py` | R1: tissue-stratified DPT sensitivity (retina / optic nerve separately) | `results/03_mg/dpt_tissue_stratified_state_medians.csv`, Fig. S4 |
| `19_gse304931_module_scores_4groups.py` | R1: per-sample module scores, four genotype-by-condition groups | `results/05_cci/GSE304931_module_scores_4groups.csv`, Fig. S5 |
| `20_spp1_receptor_benchmark.py` | R1: SPP1-receptor benchmark + per-method convergence + expression | `results/05_cci/TableS2-S4`, Fig. S6 |
| `dl_finngen.sh` | FinnGen H7_GLAUCOMA summary-statistics downloader (watchdog) | `data/` |

## Key result files cited in the manuscript

| File | Content |
| --- | --- |
| `results/02_cluster/celltype_proportion_stats.csv` | Per-sample proportion tests, 14 immune classes (retina + optic nerve) |
| `results/03_mg/mg_state_proportion_stats.csv` | Proportion tests, 12 microglial states |
| `results/03_mg/dpt_root_sensitivity.csv` | DPT root-choice sensitivity (15 roots; median ρ = 0.44, range 0.22–0.57) |
| `results/05_cci/liana_all_interactions.csv` | All candidate ligand–receptor interactions (103,619 tests) |
| `results/05_cci/liana_focus_axes.csv` | Focused SPP1–CD44 axes, OHT-specific rewiring |
| `results/06_mr_drug/mr_gene_summary.csv` | Two-sample MR summary (9 instruments / 8 genes) |
| `results/06_mr_drug/drug_reversal_ranking.csv` | LINCS L1000 reversal ranking |
| `results/03_mg/dpt_tissue_stratified_state_medians.csv` | R1: state-median DPT per tissue (Spearman ρ = 0.89, p = 8.4e-05) |
| `results/05_cci/GSE304931_module_scores_4groups.csv` | R1: per-sample module scores (6 modules x 4 groups x 4 mice) |
| `results/05_cci/GSE304931_baseline_KO_vs_WT.csv` | R1: baseline genotype tests |
| `results/05_cci/GSE304931_IOP_effect_by_genotype.csv` | R1: IOP effect within each genotype |
| `results/05_cci/TableS2_SPP1_receptor_benchmark.csv` | R1: SPP1-receptor benchmark (CD44 leading on lymphocytes) |
| `results/05_cci/TableS3_SPP1_CD44_per_method.csv` | R1: five-algorithm convergence for top SPP1-CD44 edges |
| `results/05_cci/TableS4_Spp1_Cd44_expression.csv` | R1: retina-only Spp1/Cd44 expression (mean log1p of CP10K and percent expressing) for all 14 classes by condition, with cell counts |
| `results/plot_data/` | Source data behind every figure panel |

## Changes in revision R1

- New sensitivity analysis: diffusion pseudotime recomputed independently in
  retina and optic nerve (script 18, Supplementary Fig. S4); state ordering is
  preserved (Spearman rho = 0.89, p = 8.4e-05, n = 12 states).
- New sample-level module-score analysis of GSE304931 with baseline genotype
  tests and a DAM list excluding Spp1 (script 19, Supplementary Fig. S5).
- New benchmark of annotated SPP1 receptors; CD44 carries the highest mean
  interaction magnitude among myeloid-to-lymphocyte edges (12/20 edges with
  specificity_rank <= 0.05; script 20, Supplementary Fig. S6, Tables S2-S4).
- Revised Fig. 6a/6b panels (`Fig6a_Spp1_Cd44_heatmap_revised.png`,
  `Fig6b_keygene_dotplot_revised.png`); manuscript figure numbering in R1.
- `requirements.txt` extended (leidenalg, python-igraph).

## Reproduce

```bash
pip install -r requirements.txt   # scanpy 1.12.4, liana, scrublet, ...

# 1. download GSE293358 (13 samples: barcodes/features/matrix .gz) from GEO
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293358
python3 scripts/01_qc.py
python3 scripts/02_cluster.py
python3 scripts/03_immune_recluster.py
python3 scripts/04_annotate_fig1.py
python3 scripts/05_mg_trajectory.py
python3 scripts/17_dpt_root_sensitivity.py   # root-sensitivity check

# 2. validation cohorts, communication, genetics (see table above)
```

## Data availability

- GSE293358, GSE241782, GSE306785, GSE304931 — NCBI GEO
- GTEx v8 — https://gtexportal.org; FinnGen R10 (H7_GLAUCOMA) — https://www.finngen.fi
- LINCS L1000 — via Enrichr (https://maayanlab.cloud/Enrichr)
- Processed `.h5ad` intermediates are available from the corresponding author
  upon reasonable request; they can be rebuilt from GEO raw data with
  scripts 01–05.

## Contact

Corresponding author: Xian Yang (yangxian_zhao@qdu.edu.cn).
Questions about the code: please open a GitHub issue.

## License

MIT (code). Public datasets retain their original terms of use.
