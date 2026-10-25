"""阶段6c: 双样本 MR — GTEx eQTL(暴露) vs FinnGen 青光眼(结局)
Wald ratio (单工具) / IVW (多工具); FinnGen PheWeb API 逐变异体查询
"""
import pandas as pd, numpy as np, requests, time, os, json
from scipy import stats

OUT = '/tmp/results/06_mr_drug'; os.makedirs(OUT, exist_ok=True)
iv = pd.read_csv('/tmp/mr_instruments.csv')
iv['chrom_n'] = iv.chrom.str.replace('chr', '')

rows = []
for _, r in iv.iterrows():
    url = f"https://r10.finngen.fi/api/variant/{r.chrom_n}-{r.pos}-{r.ref}-{r.alt}"
    for attempt in range(3):
        try:
            js = requests.get(url, timeout=60).json()
            break
        except Exception as e:
            print(f'retry {r.variant_id}: {e}', flush=True); time.sleep(5)
    else:
        continue
    phenos = js.get('results', []) or js.get('phenos', [])
    hit = [p for p in phenos if p.get('phenocode') == 'H7_GLAUCOMA']
    if not hit:
        print(f"{r.gene_symbol} {r.variant_id}: H7_GLAUCOMA not in {len(phenos)} phenos", flush=True)
        continue
    h = hit[0]
    rows.append({
        'gene': r.gene_symbol, 'variant_id': r.variant_id, 'chrom': r.chrom_n, 'pos': r.pos,
        'ref': r.ref, 'alt': r.alt,
        'beta_exp': r.slope, 'se_exp': r.slope_se, 'p_exp': r.pval_nominal,
        'beta_out': h.get('beta'), 'se_out': h.get('sebeta'), 'p_out': h.get('pval'),
        'maf_out': h.get('maf'), 'n_out': h.get('num_samples'),
    })
    print(f"{r.gene_symbol}: beta_exp={r.slope:.3f} beta_out={h.get('beta')} p_out={h.get('pval')}", flush=True)
    time.sleep(1)

mr = pd.DataFrame(rows)
# Wald ratio: beta_MR = beta_out / beta_exp; se = se_out / |beta_exp|
mr['beta_mr'] = mr.beta_out / mr.beta_exp
mr['se_mr'] = mr.se_out / mr.beta_exp.abs()
mr['z'] = mr.beta_mr / mr.se_mr
mr['p_mr'] = 2 * stats.norm.sf(mr.z.abs())
mr['OR'] = np.exp(mr.beta_mr)
mr['OR_lo'] = np.exp(mr.beta_mr - 1.96 * mr.se_mr)
mr['OR_hi'] = np.exp(mr.beta_mr + 1.96 * mr.se_mr)
mr.to_csv(f'{OUT}/mr_per_snp.csv', index=False)

# 按基因汇总: 单SNP=Wald; 多SNP=IVW固定效应
summ = []
for g, gdf in mr.groupby('gene'):
    if len(gdf) == 1:
        r = gdf.iloc[0]
        summ.append({'gene': g, 'n_iv': 1, 'method': 'Wald ratio',
                     'beta': r.beta_mr, 'se': r.se_mr, 'p': r.p_mr})
    else:
        w = 1 / gdf.se_mr**2
        b = (gdf.beta_mr * w).sum() / w.sum()
        se = np.sqrt(1 / w.sum())
        z = b / se
        summ.append({'gene': g, 'n_iv': len(gdf), 'method': 'IVW (FE)',
                     'beta': b, 'se': se, 'p': 2 * stats.norm.sf(abs(z))})
summ = pd.DataFrame(summ).sort_values('p')
summ['OR'] = np.exp(summ.beta)
summ['OR_lo'] = np.exp(summ.beta - 1.96*summ.se)
summ['OR_hi'] = np.exp(summ.beta + 1.96*summ.se)
summ.to_csv(f'{OUT}/mr_gene_summary.csv', index=False)
print('\n=== MR 基因汇总 (暴露=基因表达, 结局=青光眼 FinnGen H7_GLAUCOMA) ===', flush=True)
print(summ.to_string(index=False), flush=True)
print('\nDONE', flush=True)
