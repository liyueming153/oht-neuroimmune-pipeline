"""阶段5 Part2: GSE304931 Spp1-KO 小胶质 bulk
设计: WT/Spp1KO x 基线/IOP高眼压/ONC视神经挤压 (分选小胶质)
核心因果检验: DAM 轴基因在 WT+IOP 升高后, 在 Spp1KO+IOP 是否被消除?
"""
import pandas as pd, numpy as np, os, mygene
from scipy import stats

OUT = '/tmp/results/05_cci'; os.makedirs(OUT, exist_ok=True)
d = pd.read_csv('/mnt/agents/output/oht_project/data/GSE304931.csv.gz')
d['ensembl'] = d['gene'].str.split('.').str[0]
print('rows', len(d), flush=True)

# Ensembl -> symbol 映射
mg = mygene.MyGeneInfo()
ids = d['ensembl'].tolist()
res = mg.querymany(ids, scopes='ensembl.gene', fields='symbol', species='mouse',
                   as_dataframe=True, verbose=False)
mapd = res['symbol'].dropna().to_dict() if 'symbol' in res.columns else {}
d['symbol'] = d['ensembl'].map(mapd)
print('mapped:', d.symbol.notna().sum(), flush=True)
d = d.dropna(subset=['symbol'])
# 去重: 同 symbol 取均值
grp_cols = [c for c in d.columns if c not in ('gene','ensembl','symbol')]
d = d.groupby('symbol')[grp_cols].mean()
print('unique symbols:', len(d), flush=True)

GROUPS = {
    'WT': ['WT_1','WT_2','WT_3','WT_4'],
    'IOP': ['IOP_1','IOP_2','IOP_3','IOP_4'],
    'ONC': ['ONC_1','ONC_2','ONC_3'],
    'KO': ['Spp1KO_1','Spp1KO_2','Spp1KO_3','Spp1KO_4'],
    'KO_IOP': ['Spp1KO_IOP_1','Spp1KO_IOP_2','Spp1KO_IOP_3','Spp1KO_IOP_4'],
    'KO_ONC': ['Spp1KO_ONC_1','Spp1KO_ONC_2','Spp1KO_ONC_3','Spp1KO_ONC_4'],
}
# 值为连续浮点, 视为 normalized 计数; 转 log2(x+1) 空间
L = np.log2(d[ [c for g in GROUPS.values() for c in g] ] + 1)

def de(g1, g2, name):
    a = L[GROUPS[g1]].values; b = L[GROUPS[g2]].values
    fc = a.mean(1) - b.mean(1)          # log2 空间差 = log2FC
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy='omit')
    out = pd.DataFrame({'gene': L.index, 'log2FC': fc, 'pvalue': p}).dropna()
    out = out.sort_values('pvalue')
    out['padj'] = stats.false_discovery_control(out['pvalue']) if hasattr(stats,'false_discovery_control') else np.minimum.accumulate(out['pvalue'].values*len(out)/np.arange(1,len(out)+1))
    out.to_csv(f'{OUT}/GSE304931_DE_{name}.csv', index=False)
    print(f'{name}: {len(out)} genes, FDR<0.05: {(out.padj<0.05).sum()}, '
          f'FDR<0.05&|FC|>1: {((out.padj<0.05)&(out.log2FC.abs()>1)).sum()}', flush=True)
    return out.set_index('gene')

de_iop = de('IOP', 'WT', 'IOP_vs_WT')
de_koiop = de('KO_IOP', 'KO', 'KO_IOP_vs_KO')
de_ko_vs_iop = de('KO_IOP', 'IOP', 'KO_IOP_vs_WT_IOP')
de_onc = de('ONC', 'WT', 'ONC_vs_WT')

# ===== 阶段2/3 DAM/IFN 轴签名基因在 Spp1-KO 下的依赖性 =====
SIG = {
 'DAM_MG': ['Apoe','Spp1','Lgals3','Cst7','Lpl','Ctsb','Ctsd','Ctsl','Trem2','Axl','Itgax','Cd9','Lilrb4','Clec7a','Gpnmb','Fabp5','Lgals1','Csf1','B2m','Fth1'],
 'IFN_module': ['Ifit3','Isg15','Usp18','Irf7','Ifit1','Ifi44','Oasl2','Rsad2','Stat1','Mx1','Ifit2','Oas1a','Bst2','Ifi27l2a'],
 'MHCII': ['Cd74','H2-Aa','H2-Ab1','H2-Eb1','H2-DMa','H2-DMb1'],
 'Homeostatic_MG': ['P2ry12','Tmem119','Sall1','Hexb','Selplg','Cx3cr1','Siglech','Mafb','Csf1r','Tgfbr1'],
 'Chemokine': ['Ccl2','Ccl3','Ccl4','Ccl5','Ccl12','Cxcl10','Cxcl2','Cxcl16'],
}
# 每个签名: 方向一致性 + WT-IOP变化 vs KO-IOP变化 的配对比较
rows = []
for sig, genes in SIG.items():
    genes = [g for g in genes if g in L.index]
    r = {'signature': sig, 'n_genes': len(genes)}
    for nm, tbl in [('WT_IOP', de_iop), ('KO_IOP', de_koiop), ('WT_ONC', de_onc)]:
        sub = tbl.reindex(genes)
        r[f'{nm}_meanFC'] = sub.log2FC.mean()
        r[f'{nm}_frac_up'] = (sub.log2FC > 0).mean()
        r[f'{nm}_paired_p'] = stats.wilcoxon(sub.log2FC).pvalue if len(sub.dropna())>=6 else np.nan
    # Spp1 依赖性: WT-IOP 的 FC 是否显著大于 KO-IOP 的 FC
    both = pd.concat([de_iop.reindex(genes).log2FC.rename('wt'),
                      de_koiop.reindex(genes).log2FC.rename('ko')], axis=1).dropna()
    r['Spp1_dependence_p'] = stats.wilcoxon(both.wt, both.ko).pvalue if len(both)>=6 else np.nan
    r['Spp1_dependence_dFC'] = (both.wt - both.ko).mean()
    rows.append(r)
sig_df = pd.DataFrame(rows)
sig_df.to_csv(f'{OUT}/GSE304931_signature_dependence.csv', index=False)
print('\n=== 签名依赖性 ===', flush=True)
print(sig_df.to_string(index=False), flush=True)

# 关键基因明细
key = ['Spp1','Apoe','Lgals3','Cst7','Itgax','Gpnmb','Trem2','Ifit3','Isg15','Usp18',
       'Cd74','H2-Aa','Ccl5','Cxcl10','C3','Cx3cl1','Cx3cr1','P2ry12','Tmem119']
kr = []
for g in key:
    if g not in L.index: continue
    row = {'gene': g,
           'mean_WT': L.loc[g, GROUPS['WT']].mean(), 'mean_IOP': L.loc[g, GROUPS['IOP']].mean(),
           'mean_KO': L.loc[g, GROUPS['KO']].mean(), 'mean_KO_IOP': L.loc[g, GROUPS['KO_IOP']].mean()}
    for nm, tbl in [('IOP_vs_WT', de_iop), ('KO_IOP_vs_KO', de_koiop), ('KO_IOP_vs_WT_IOP', de_ko_vs_iop)]:
        if g in tbl.index:
            row[f'{nm}_FC'] = tbl.loc[g, 'log2FC']; row[f'{nm}_p'] = tbl.loc[g, 'pvalue']
    kr.append(row)
pd.DataFrame(kr).to_csv(f'{OUT}/GSE304931_keygenes.csv', index=False)
print('\n=== 关键基因 ===', flush=True)
print(pd.DataFrame(kr).to_string(index=False), flush=True)
print('\nDONE', flush=True)
