#!/usr/bin/env python3
"""Step 09: Mendelian randomization (Methods 4.7).

Instruments: significant GTEx v8 whole-blood cis-eQTLs for eight core genes
(SPP1, CX3CR1, TREM2, CD44, ITGAX, CCL5, CSF1R, LGALS3), pruned at 250 kb
(minimum p per window, <=10 per gene; p<5e-8 preferred).
Outcome: FinnGen R10 H7_GLAUCOMA (20,906 cases, 391,275 controls).
Harmonization on the alternative allele; Wald ratio for single instruments,
fixed-effect IVW for multi-instrument genes. OR per 1-SD expression.

Requires: gtex_eqtl.csv (rsid, gene, beta_exp, se_exp) and
finngen_h7_glaucoma.csv (rsid, beta_out, se_out) placed in
/mnt/agents/work/mr/.
"""
import numpy as np
import os
import pandas as pd
from scipy import stats

DATA = os.environ.get('MR_DIR', 'data/MR')
OUT = os.environ.get('OUT_DIR', 'out')

GENES = ['SPP1', 'CX3CR1', 'TREM2', 'CD44', 'ITGAX', 'CCL5', 'CSF1R', 'LGALS3']
exp = pd.read_csv(f'{DATA}/gtex_eqtl.csv')
out = pd.read_csv(f'{DATA}/finngen_h7_glaucoma.csv')

exp = exp[exp.gene.isin(GENES)]
exp = exp.sort_values('p_exp' if 'p_exp' in exp else 'se_exp')
exp = exp.groupby('gene').head(10)  # <=10 per gene after 250-kb pruning
harm = exp.merge(out, on='rsid')

rows = []
for gene, g in harm.groupby('gene'):
    if len(g) == 1:
        b = g.beta_out.iloc[0] / g.beta_exp.iloc[0]
        se = abs(g.se_out.iloc[0] / g.beta_exp.iloc[0])
        method = 'Wald ratio'
    else:
        w = 1 / g.se_out**2
        b = np.sum(w * g.beta_out / g.beta_exp) / np.sum(w)
        se = np.sqrt(1 / np.sum(w))
        method = 'IVW (fixed effect)'
    z = b / se
    rows.append(dict(gene=gene, n_instruments=len(g), method=method,
                     OR=np.exp(b), OR_lo=np.exp(b - 1.96 * se),
                     OR_hi=np.exp(b + 1.96 * se), p=2 * stats.norm.sf(abs(z))))
res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/mr_results.csv', index=False)
print(res.to_string(index=False))
