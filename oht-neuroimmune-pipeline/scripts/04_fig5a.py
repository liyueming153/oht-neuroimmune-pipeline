#!/usr/bin/env python3
"""Step 04: Regenerate Figure 5a - Spp1->CD44 sender x receiver heatmaps
(Sham vs MB, Retina) from the 1000-permutation LIANA results.
Cells failing the specificity criterion are masked (shown white)."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = os.environ.get('OUT_DIR', 'out')
NPERMS = 1000
SPEC_CUT = 0.05  # specificity_rank criterion (LIANA default for rank_aggregate)

sham = pd.read_csv(f'{OUT}/liana_Retina_Sham_{NPERMS}perms.csv.gz')
mb = pd.read_csv(f'{OUT}/liana_Retina_MB_{NPERMS}perms.csv.gz')

ORDER = ['Homeostatic MG', 'Activated MG (Ccl5+)', 'DAM-like MG',
         'Iron+ MG (Fth1+)', 'Proliferating MG', 'BAM', 'MHCII+ Mac/APC',
         'pDC', 'T cell', 'gdT cell', 'NK cell', 'B cell', 'Neutrophil',
         'Monocyte']

def mat(df):
    d = df[(df.ligand_complex == 'Spp1') & (df.receptor_complex == 'Cd44')].copy()
    d.loc[d.specificity_rank > SPEC_CUT, 'lr_means'] = np.nan
    labels = [x for x in ORDER if x in set(d.source) | set(d.target)
              or x in set(df.source) | set(df.target)]
    M = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=float)
    for _, r in d.iterrows():
        M.loc[r.source, r.target] = r.lr_means
    return M, labels

M_sham, labels = mat(sham)
M_mb, _ = mat(mb)
# unify label sets
labels = [x for x in ORDER if x in set(M_sham.index) | set(M_mb.index)]
M_sham = M_sham.reindex(index=labels, columns=labels)
M_mb = M_mb.reindex(index=labels, columns=labels)

fig, axes = plt.subplots(1, 2, figsize=(19.84, 9.0), dpi=150)
vmax = np.nanmax([np.nanmax(M_sham.values), np.nanmax(M_mb.values), 1.0])
for ax, M, title in zip(axes, [M_sham, M_mb],
                        ['Spp1→CD44 (Sham, Retina)', 'Spp1→CD44 (MB, Retina)']):
    im = ax.imshow(M.values, cmap='Reds', vmin=0, vmax=vmax,
                   aspect='equal')
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=90, fontsize=13)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=13)
    ax.set_title(f'{title}\nsender=row, receiver=col', fontsize=16)
    ax.set_xlabel('receiver', fontsize=15)
    if ax is axes[0]:
        ax.set_ylabel('sender', fontsize=15)
fig.subplots_adjust(left=0.09, right=0.90, top=0.86, bottom=0.22, wspace=0.28)
cax = fig.add_axes([0.915, 0.30, 0.012, 0.45])
cb = fig.colorbar(im, cax=cax)
cb.ax.tick_params(labelsize=12)
cb.set_label('LR means', fontsize=14)
fig.savefig(f'{OUT}/fig5a_1000perms.png', dpi=150)
print('saved', f'{OUT}/fig5a_1000perms.png')
print('Sham significant pairs:', int((M_sham.notna()).sum().sum()))
print('MB significant pairs:', int((M_mb.notna()).sum().sum()))
top = M_mb.stack().sort_values(ascending=False).head(5)
print(top)
