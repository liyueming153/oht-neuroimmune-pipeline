"""阶段6b: 从 GTEx 全血 signif pairs 提取核心基因 cis-eQTL 工具变量"""
import pandas as pd, numpy as np, subprocess, io

raw = subprocess.run(['zcat', '/tmp/wb_signif_pairs.txt.gz'], capture_output=True, text=True).stdout
df = pd.read_csv(io.StringIO(raw), sep='\t')
print('total pairs:', len(df), 'cols:', df.columns.tolist()[:12], flush=True)

GENES = {
 'SPP1':'ENSG00000118785','CX3CR1':'ENSG00000168329','CX3CL1':'ENSG00000006210',
 'CD44':'ENSG00000026508','TREM2':'ENSG00000095970','LGALS3':'ENSG00000131981',
 'CCL5':'ENSG00000271503','ITGAX':'ENSG00000140678','APOE':'ENSG00000130203',
 'CSF1R':'ENSG00000182578','TGFB1':'ENSG00000105329','IL34':'ENSG00000157368',
 'C3':'ENSG00000125730','CCL2':'ENSG00000108691','CXCL10':'ENSG00000169245',
}
df['ens'] = df['gene_id'].str.split('.').str[0]
df = df[df.ens.isin(GENES.values())]
print('pairs in target genes:', len(df), flush=True)

# 解析 variant_id: chr1_12345_A_G_b38
v = df['variant_id'].str.split('_', expand=True)
df['chrom'] = v[0]; df['pos'] = v[1].astype(int); df['ref'] = v[2]; df['alt'] = v[3]

# 工具变量选择: p<5e-8 优先, 不足则放宽到 1e-6; 250kb 窗口去连锁(保留最小p)
def select(gdf):
    gdf = gdf.sort_values('pval_nominal')
    for thr in [5e-8, 1e-6]:
        cand = gdf[gdf.pval_nominal < thr]
        if len(cand) >= 2: break
    cand = gdf if len(cand) < 2 else cand
    kept = []
    used = cand.sort_values('pval_nominal')
    for _, r in used.iterrows():
        if all(abs(r.pos - k) > 250000 for k in kept):
            kept.append(r.pos)
        if len(kept) >= 10: break
    return used[used.pos.isin(kept)]

iv = df.groupby('ens', group_keys=False).apply(select)
iv['gene_symbol'] = iv.ens.map({v:k for k,v in GENES.items()})
iv = iv[['gene_symbol','variant_id','chrom','pos','ref','alt','slope','slope_se','pval_nominal','tss_distance','maf']]
iv.to_csv('/tmp/mr_instruments.csv', index=False)
print('\n=== 每个基因的工具变量数 ===', flush=True)
print(iv.groupby('gene_symbol').size().sort_values(ascending=False).to_string(), flush=True)
print('\ntotal IVs:', len(iv), flush=True)
print(iv[['gene_symbol','variant_id','slope','slope_se','pval_nominal']].head(30).to_string(index=False), flush=True)
