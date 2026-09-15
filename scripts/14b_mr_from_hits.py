"""阶段6c: MR 计算 — GTEx eQTL(暴露) + FinnGen H7_GLAUCOMA(结局, 本地hits)
两侧均按 alt 等位基因对齐; Wald ratio / IVW
"""
import pandas as pd, numpy as np, os
from scipy import stats

OUT = '/tmp/results/06_mr_drug'; os.makedirs(OUT, exist_ok=True)
iv = pd.read_csv('/tmp/mr_instruments.csv')
iv['chrom_n'] = iv.chrom.str.replace('chr','').astype(str)
iv['key'] = iv.chrom_n + ':' + iv.pos.astype(str)

hits = pd.read_csv('/tmp/finngen_hits.tsv', sep='\t',
                   names=['chrom','pos','ref','alt','rsids','nearest','pval','mlogp','beta','sebeta','af_alt','af_case','af_ctrl'])
hits['chrom'] = hits['chrom'].astype(str)
hits['key'] = hits.chrom + ':' + hits.pos.astype(str)

m = iv.merge(hits, on='key', suffixes=('_exp','_out'))
# 等位基因核对
m['allele_match'] = (m.ref_exp == m.ref_out) & (m.alt_exp == m.alt_out)
print('allele match:', m.allele_match.sum(), '/', len(m), flush=True)
m = m[m.allele_match]

m['beta_mr'] = m.beta / m.slope
m['se_mr'] = m.sebeta / m.slope.abs()
m['z'] = m.beta_mr / m.se_mr
m['p_mr'] = 2 * stats.norm.sf(m.z.abs())
m['OR'] = np.exp(m.beta_mr)
m['OR_lo'] = np.exp(m.beta_mr - 1.96*m.se_mr)
m['OR_hi'] = np.exp(m.beta_mr + 1.96*m.se_mr)
out = m[['gene_symbol','variant_id','rsids','slope','slope_se','pval_nominal','beta','sebeta','pval',
         'beta_mr','se_mr','p_mr','OR','OR_lo','OR_hi']]
out.columns = ['gene','variant','rsid','beta_exp','se_exp','p_exp','beta_out','se_out','p_out',
               'beta_mr','se_mr','p_mr','OR','OR_lo','OR_hi']
out.to_csv(f'{OUT}/mr_per_snp.csv', index=False)

summ = []
for g, gdf in m.groupby('gene_symbol'):
    if len(gdf) == 1:
        r = gdf.iloc[0]
        summ.append({'gene': g, 'n_iv': 1, 'method': 'Wald ratio', 'beta': r.beta_mr,
                     'se': r.se_mr, 'p': r.p_mr, 'iv_strength_p': r.pval_nominal})
    else:
        w = 1/gdf.se_mr**2
        b = (gdf.beta_mr*w).sum()/w.sum(); se = np.sqrt(1/w.sum())
        summ.append({'gene': g, 'n_iv': len(gdf), 'method': 'IVW (FE)', 'beta': b, 'se': se,
                     'p': 2*stats.norm.sf(abs(b/se)), 'iv_strength_p': gdf.pval_nominal.min()})
summ = pd.DataFrame(summ).sort_values('p')
summ['OR'] = np.exp(summ.beta)
summ['OR_lo'] = np.exp(summ.beta - 1.96*summ.se)
summ['OR_hi'] = np.exp(summ.beta + 1.96*summ.se)
summ['FDR'] = stats.false_discovery_control(summ.p)
summ.to_csv(f'{OUT}/mr_gene_summary.csv', index=False)
cols = ['gene','n_iv','method','beta','se','p','FDR','OR','OR_lo','OR_hi']
print('\n=== MR 汇总: 基因表达 → 青光眼 (FinnGen R10 H7_GLAUCOMA, 20906 cases / 391275 ctrls) ===', flush=True)
print(summ[cols].to_string(index=False, float_format=lambda x: f'{x:.4g}'), flush=True)
print('\n逐SNP明细:', flush=True)
print(out.to_string(index=False, float_format=lambda x: f'{x:.4g}'), flush=True)
print('\nDONE', flush=True)
