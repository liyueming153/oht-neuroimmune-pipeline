"""阶段5 Part1: 细胞通讯分析 (liana, mouseconsensus 资源)
视网膜 Sham vs MB 分开跑, 聚焦 Cx3cl1-Cx3cr1 / Spp1 轴 / 补体 / 趋化因子
"""
import scanpy as sc, liana as li, pandas as pd, numpy as np, gc, os, warnings
warnings.filterwarnings('ignore')

OUT = '/tmp/results/05_cci'; os.makedirs(OUT, exist_ok=True)
print('loading h5ad...', flush=True)
a = sc.read_h5ad('/tmp/GSE293358_final.h5ad')
print('loaded', a.shape, flush=True)

# liana 需要 log-normalized 的 X (当前即为 log1p norm 状态)
res_all = {}
for tissue in ['Retina', 'OpticNerve']:
    for cond in ['Sham', 'MB']:
        sub = a[(a.obs.tissue == tissue) & (a.obs.condition == cond)]
        print(f'=== {tissue} {cond}: {sub.n_obs} cells ===', flush=True)
        if sub.n_obs < 2000:
            print('skip (too few)', flush=True); continue
        sub = sub.copy()
        gc.collect()
        li.mt.rank_aggregate(
            adata=sub, groupby='celltype',
            resource_name='mouseconsensus',
            expr_prop=0.1, min_cells=10, n_perms=100,
            use_raw=False, verbose=False, inplace=True,
            key_added='liana_res')
        df = sub.uns['liana_res'].copy()
        df['tissue'] = tissue; df['condition'] = cond
        res_all[f'{tissue}_{cond}'] = df
        del sub; gc.collect()
        print(f'  -> {len(df)} interactions', flush=True)

res = pd.concat(res_all, ignore_index=True)
res.to_csv(f'{OUT}/liana_all_interactions.csv', index=False)

# 聚焦关键轴
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
    m = res[(res.ligand == lg) & (res.receptor == rc)]
    for _, r in m.iterrows():
        rows.append({'axis': name, 'tissue': r.tissue, 'condition': r.condition,
                     'source': r.source, 'target': r.target,
                     'lr_means': r.lr_means, 'magnitude_rank': r.magnitude_rank,
                     'specificity_rank': r.specificity_rank})
focus = pd.DataFrame(rows)
focus.to_csv(f'{OUT}/liana_focus_axes.csv', index=False)
print('\n聚焦轴互作数:', len(focus), flush=True)

# Cx3cl1-Cx3cr1 全景: 谁在发 (source), 谁在收 (target)
cx = focus[focus.axis == 'Cx3cl1|Cx3cr1'].copy()
print('\n=== Cx3cl1-Cx3cr1 (Retina) ===', flush=True)
print(cx[cx.tissue == 'Retina'].sort_values(['condition', 'lr_means'], ascending=[True, False])
      [['condition','source','target','lr_means','magnitude_rank','specificity_rank']].to_string(index=False), flush=True)

sp = focus[focus.axis.str.startswith('Spp1')].copy()
print('\n=== Spp1 轴 (Retina) ===', flush=True)
print(sp[sp.tissue == 'Retina'].sort_values(['axis','condition','lr_means'], ascending=[True,True,False])
      [['axis','condition','source','target','lr_means','magnitude_rank','specificity_rank']].head(60).to_string(index=False), flush=True)

# 配体/受体基因在各细胞类型的表达 (用于画图)
genes = sorted({g for pair in AXES.values() for g in pair})
print('\n关键基因在 var_names 中存在:', [g for g in genes if g in a.var_names], flush=True)
expr = []
for tissue in ['Retina']:
    for ct in a.obs.celltype.cat.categories:
        for cond in ['Sham', 'MB']:
            m = (a.obs.tissue == tissue) & (a.obs.celltype == ct) & (a.obs.condition == cond)
            if m.sum() < 20: continue
            idx = np.where(np.asarray(m))[0]
            row = {'tissue': tissue, 'celltype': ct, 'condition': cond, 'n': int(m.sum())}
            for g in genes:
                if g in a.var_names:
                    j = a.var_names.get_loc(g)
                    # 分块取避免列切片CSC问题: 取单列用 .X[:, j] 稀疏列 -> dense
                    v = a.X[idx][:, j]
                    row[g] = float(np.asarray(v.mean(axis=0)).ravel()[0])
                    row[g + '_pct'] = float((v > 0).mean())
            expr.append(row)
pd.DataFrame(expr).to_csv(f'{OUT}/key_gene_expr_by_celltype.csv', index=False)
print('\nDONE', flush=True)
