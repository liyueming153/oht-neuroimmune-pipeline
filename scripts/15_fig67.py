"""阶段6: Fig6 MR 森林图 + Fig7 药物重定位图"""
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/tmp/results/06_mr_drug'
plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.8})

# ---- Fig6: MR 森林图 ----
s = pd.read_csv(f'{OUT}/mr_gene_summary.csv').sort_values('beta')
fig, ax = plt.subplots(figsize=(6.2, 4.2))
y = np.arange(len(s))
ax.errorbar(s.OR, y, xerr=[s.OR-s.OR_lo, s.OR_hi-s.OR], fmt='o', color='#4C78A8',
            ecolor='#4C78A8', elinewidth=1.5, capsize=3, markersize=6)
ax.axvline(1, color='k', ls='--', lw=0.8)
for i, (_, r) in enumerate(s.iterrows()):
    star = '*' if r.p < 0.05 else ''
    ax.text(max(s.OR_hi)*1.02, i, f'OR={r.OR:.2f} ({r.OR_lo:.2f}-{r.OR_hi:.2f}), p={r.p:.3f}{star}',
            va='center', fontsize=7.5)
ax.set_yticks(y, [f'{g} (n_iv={n})' for g, n in zip(s.gene, s.n_iv)], fontsize=8.5)
ax.set_xlabel('OR for glaucoma per 1-SD higher blood expression (95% CI)')
ax.set_xlim(min(s.OR_lo)*0.92, max(s.OR_hi)*1.9)
ax.set_title('Two-sample MR: GTEx whole-blood cis-eQTL → FinnGen R10 glaucoma\n(20,906 cases / 391,275 controls)', fontsize=9)
plt.tight_layout(); plt.savefig(f'{OUT}/Fig6_MR_forest.png', dpi=300, bbox_inches='tight'); plt.close()
print('Fig6 done', flush=True)

# ---- Fig7: 药物重定位 ----
rev = pd.read_csv(f'{OUT}/drug_reversal_ranking.csv').head(15).iloc[::-1]
mech = {
 'AS-601245':'JNK inhibitor','TPCA-1':'IKK2/NF-kB inhibitor','BI-2536':'PLK1 inhibitor',
 'AZD-8055':'mTOR inhibitor','BMS-387032':'CDK2/7/9 inhibitor','radicicol':'HSP90 inhibitor',
 'crizotinib':'c-MET/ALK inhibitor (FDA-approved)','SB590885':'BRAF inhibitor',
 'geldanamycin':'HSP90 inhibitor','buparlisib':'pan-PI3K inhibitor','withaferin-a':'steroidal lactone (NF-kB)',
 'tivantinib':'c-MET inhibitor','CGP-60474':'CDK inhibitor','QL-XII-47':'BTK inhibitor','PI-103':'PI3K/mTOR inhibitor'}
fig, ax = plt.subplots(figsize=(7.2, 4.8))
colors = plt.cm.YlOrRd(np.linspace(0.35, 0.85, len(rev)))
ax.barh(rev.compound, rev.reversal_score, color=colors)
for i, (_, r) in enumerate(rev.iterrows()):
    ax.text(r.reversal_score + 0.15, i, mech.get(r.compound, ''), va='center', fontsize=7.2, color='#555555')
ax.set_xlabel('Signature reversal score (-log10 adjP, weighted across 4 modules)')
ax.set_title('Drug repurposing: LINCS L1000 compounds reversing the OHT immune signature\n(DAM/IFN/Chemokine up-modules in DOWN library; Homeostatic module in UP library)', fontsize=9)
ax.set_xlim(0, rev.reversal_score.max()*1.55)
plt.tight_layout(); plt.savefig(f'{OUT}/Fig7_drug_reversal.png', dpi=300, bbox_inches='tight'); plt.close()
print('Fig7 done', flush=True)
print('ALL DONE', flush=True)
