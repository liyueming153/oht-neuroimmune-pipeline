#!/usr/bin/env python3
"""Step 07: Sorted microglia Spp1-KO cohort (GSE304931; Methods 4.5).

Ensembl IDs mapped to symbols (MyGene); analyzed as log2(normalized+1);
per-gene contrasts use Welch t-tests between independent animal groups;
module-level summaries compare per-gene fold-changes of module genes between
the two contrasts (paired Wilcoxon across genes; one-sided Mann-Whitney U for
the dependence contrast) -- gene-level descriptive summaries, not animal-level
inference (see Methods).
"""
import numpy as np
import os
import pandas as pd
from scipy import stats

DATA = os.environ.get('GSE304931_DIR', 'data/GSE304931')
OUT = os.environ.get('OUT_DIR', 'out')

df = pd.read_csv(f'{DATA}/GSE304931_normalized_counts.csv.gz', index_col=0)
# columns follow the GEO deposit: WT, WT+IOP, KO, KO+IOP, (+ONC arms)
groups = {  # column lists per group; adjust to deposited column names
    'WT':      [c for c in df.columns if c.startswith('WT_')],
    'WT_IOP':  [c for c in df.columns if c.startswith('WTIOP_')],
    'KO':      [c for c in df.columns if c.startswith('KO_')],
    'KO_IOP':  [c for c in df.columns if c.startswith('KOIOP_')],
}
log = np.log2(df + 1)

def contrast(a_cols, b_cols):
    rows = []
    for gene in log.index:
        a, b = log.loc[gene, a_cols].values, log.loc[gene, b_cols].values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append(dict(gene=gene, log2FC=a.mean() - b.mean(), p=p))
    return pd.DataFrame(rows)

wt_iop = contrast(groups['WT_IOP'], groups['WT'])
ko_iop = contrast(groups['KO_IOP'], groups['KO'])
ko_wt  = contrast(groups['KO_IOP'], groups['WT_IOP'])
wt_iop.to_csv(f'{OUT}/gse304931_WTIOP_vs_WT.csv', index=False)
ko_iop.to_csv(f'{OUT}/gse304931_KOIOP_vs_KO.csv', index=False)
ko_wt.to_csv(f'{OUT}/gse304931_KOIOP_vs_WTIOP.csv', index=False)

for g in ['Spp1', 'Cx3cr1', 'Cx3cl1']:
    for name, d in [('WT+IOP vs WT', wt_iop), ('KO+IOP vs WT+IOP', ko_wt)]:
        m = d.gene.str.upper() == g.upper()
        print(g, name, d.loc[m, ['log2FC', 'p']].to_string(index=False, header=False))

DAM = ['Apoe','Spp1','Lgals3','Cst7','Trem2','Tyrobp','Gpnmb','Itgax','Cd9','Ctsb','Ctsd','Csf1']
IFN = ['Ifit3','Isg15','Usp18','Ifit1','Oasl2','Irf7','Stat1']

def module_fc(d, genes):
    m = d.gene.str.upper().isin([g.upper() for g in genes])
    return d.loc[m, 'log2FC']

dam_wt = module_fc(wt_iop, DAM); dam_ko = module_fc(ko_iop, DAM)
w, p = stats.wilcoxon(dam_ko.values, dam_wt.values)
print(f'DAM dependence: KO+IOP mean {dam_ko.mean():.2f} vs WT+IOP {dam_wt.mean():.2f}, Wilcoxon p={p:.4g}')

ifn_wt = module_fc(wt_iop, IFN); ifn_ko = module_fc(ko_iop, IFN)
w, p = stats.wilcoxon(ifn_ko.values, ifn_wt.values)
print(f'IFN module: KO+IOP mean {ifn_ko.mean():.2f} vs WT+IOP {ifn_wt.mean():.2f}, Wilcoxon p={p:.4g}')
