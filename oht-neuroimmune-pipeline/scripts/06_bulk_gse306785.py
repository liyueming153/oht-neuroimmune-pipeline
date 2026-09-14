#!/usr/bin/env python3
"""Step 06: Sorted Muller glia DE and astrocyte-reactivity modules
(GSE306785; Methods 4.5). Normalized counts analyzed in log2 space
(Welch t, BH). Generates the Figure 4b module barplot values and the
Cx3cl1 / Spp1 / Gfap gene-level statistics (Figure 4c)."""
import numpy as np
import os
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

DATA = os.environ.get('GSE306785_DIR', 'data/GSE306785')
OUT = os.environ.get('OUT_DIR', 'out')

df = pd.read_csv(f'{DATA}/GSE306785_macs_normalized_counts.csv.gz', index_col=0)
sham_cols = [c for c in df.columns if 'SHAM' in c.upper()]
oht_cols = [c for c in df.columns if c not in sham_cols]
log = np.log2(df + 1)

# gene-level stats
rows = []
for gene, rec in df.iterrows():
    a, b = log.loc[gene, oht_cols].values, log.loc[gene, sham_cols].values
    t, p = stats.ttest_ind(a, b, equal_var=False)
    rows.append(dict(gene=gene, log2FC=a.mean() - b.mean(), p=p))
de = pd.DataFrame(rows)
de['fdr'] = stats.false_discovery_control(de['p'])
de.to_csv(f'{OUT}/gse306785_de.csv', index=False)
for g in ['Cx3cl1', 'Spp1', 'Gfap']:
    m = de.gene.str.upper() == g.upper()
    print(de.loc[m].to_string(index=False))

MODULES = {
 'Pan_reactive': ['Gfap','Vim','Serpina3n','Timp1','Osmr','Cd44'],
 'A1_neurotoxic': ['H2-D1','Gbp2','Serping1','C3','Fbln5','Amigo2','Psmb8'],
 'A2_neuroprotective': ['S100a10','Tgm1','Ptx3','Emp1','Clcf1'],
 'Complement': ['C1qa','C1qb','C1qc','C3','C4b'],
 'Chemokine_recruit': ['Ccl2','Ccl5','Cxcl10','Ccl12'],
}
idx = {g.upper(): i for i, g in enumerate(log.index)}
means_s, sds_s, means_o, sds_o, names = [], [], [], [], []
for name, genes in MODULES.items():
    rows_ = [idx[g.upper()] for g in genes if g.upper() in idx]
    sc_ = log.iloc[rows_].mean(axis=0)
    s, o = sc_[sham_cols].values, sc_[oht_cols].values
    names.append(name)
    means_s.append(s.mean()); sds_s.append(s.std(ddof=1))
    means_o.append(o.mean()); sds_o.append(o.std(ddof=1))
    print(f'{name}: SHAM={s.mean():.2f}±{s.std(ddof=1):.2f}  OHT={o.mean():.2f}±{o.std(ddof=1):.2f}')

x = np.arange(len(names)); w = 0.38
fig, ax = plt.subplots(figsize=(15.79, 8.8), dpi=150)
ax.bar(x - w/2, means_s, w, yerr=sds_s, capsize=8, color='#9AB0A9', label='SHAM',
       error_kw=dict(ecolor='black', lw=2.5))
ax.bar(x + w/2, means_o, w, yerr=sds_o, capsize=8, color='#C66444', label='OHT',
       error_kw=dict(ecolor='black', lw=2.5))
ax.set_xticks(x); ax.set_xticklabels(names, rotation=20, ha='right', fontsize=17)
ax.set_ylabel('Mean log2(norm+1)', fontsize=18); ax.set_ylim(0, 18)
ax.tick_params(axis='y', labelsize=15)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.8); ax.spines['bottom'].set_linewidth(1.8)
ax.yaxis.grid(True, alpha=0.3); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=17, loc='upper right')
fig.tight_layout()
fig.savefig(f'{OUT}/fig4b_modules.png', dpi=150)
print('figure saved')
