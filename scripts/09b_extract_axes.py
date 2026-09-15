"""阶段5 Part1b: 从已保存的 liana 结果提取聚焦轴 + 关键基因表达"""
import scanpy as sc, pandas as pd, numpy as np, os, warnings
warnings.filterwarnings('ignore')
OUT = '/tmp/results/05_cci'

res = pd.read_csv(f'{OUT}/liana_all_interactions.csv')
print('interactions:', len(res), res.tissue.unique(), res.condition.unique(), flush=True)

AXES = {
    'Cx3cl1|Cx3cr1': ('Cx3cl1', 'Cx3cr1'),
    'Spp1|Cd44':     ('Spp1', 'Cd44'),
    'Spp1|Itgav':    ('Spp1', 'Itgav'),
    'Spp1|Itgb1':    ('Spp1', 'Itgb1'),
    'Spp1|Itgb5':    ('Spp1', 'Itgb5'),
    'Spp1|Itga9':    ('Spp1', 'Itga9'),
    'C3|C3ar1':      ('C3', 'C3ar1'),
    'Ccl5|Ccr5':     ('Ccl5', 'Ccr5'),
    'Cxcl10|Cxcr3':  ('Cxcl10', 'Cxcr3'),
    'Il34|Csf1r':    ('Il34', 'Csf1r'),
    'Tgfb1|Tgfbr1':  ('Tgfb1', 'Tgfbr1'),
    'Gas6|Axl':      ('Gas6', 'Axl'),
    'Apoe|Trem2':    ('Apoe', 'Trem2'),
    'Apoe|Lrp1':     ('Apoe', 'Lrp1'),
}
rows = []
for name, (lg, rc) in AXES.items():
    m = res[(res.ligand_complex == lg) & (res.receptor_complex == rc)]
    for _, r in m.iterrows():
        rows.append({'axis': name, 'tissue': r.tissue, 'condition': r.condition,
                     'source': r.source, 'target': r.target,
                     'lr_means': r.lr_means, 'magnitude_rank': r.magnitude_rank,
                     'specificity_rank': r.specificity_rank})
focus = pd.DataFrame(rows)
focus.to_csv(f'{OUT}/liana_focus_axes.csv', index=False)
print('focus rows:', len(focus), flush=True)

print('\n=== Cx3cl1-Cx3cr1 (Retina) ===', flush=True)
cx = focus[(focus.axis=='Cx3cl1|Cx3cr1') & (focus.tissue=='Retina')]
print(cx.sort_values(['condition','lr_means'], ascending=[True,False])
      [['condition','source','target','lr_means','magnitude_rank']].to_string(index=False), flush=True)

print('\n=== Spp1 轴 top (Retina) ===', flush=True)
sp = focus[focus.axis.str.startswith('Spp1') & (focus.tissue=='Retina')]
print(sp.sort_values(['axis','condition','lr_means'], ascending=[True,True,False])
      [['axis','condition','source','target','lr_means','magnitude_rank']].head(50).to_string(index=False), flush=True)

# 关键基因表达 (Retina, 逐细胞类型)
a = sc.read_h5ad('/tmp/GSE293358_final.h5ad')
genes = sorted({g for p in AXES.values() for g in p})
print('\n存在于矩阵:', [g for g in genes if g in a.var_names], flush=True)
expr = []
for ct in a.obs.celltype.cat.categories:
    for cond in ['Sham', 'MB']:
        m = np.asarray((a.obs.tissue=='Retina') & (a.obs.celltype==ct) & (a.obs.condition==cond))
        if m.sum() < 20: continue
        idx = np.where(m)[0]
        row = {'celltype': ct, 'condition': cond, 'n': int(m.sum())}
        for g in genes:
            if g not in a.var_names: continue
            j = a.var_names.get_loc(g)
            v = a.X[idx][:, j]
            row[g] = float(np.asarray(v.mean(axis=0)).ravel()[0])
            row[g+'_pct'] = float(np.asarray((v>0).mean(axis=0)).ravel()[0])
        expr.append(row)
pd.DataFrame(expr).to_csv(f'{OUT}/key_gene_expr_by_celltype.csv', index=False)
print('\nDONE', flush=True)
