"""阶段6d: 药物重定位 — DAM/IFN 签名 vs Enrichr LINCS L1000
逻辑: 疾病中上调的签名基因, 若富集于某化合物的 DOWN 扰动库, 则该化合物可逆转疾病签名
"""
import requests, json, pandas as pd, numpy as np, os, time
OUT = '/tmp/results/06_mr_drug'; os.makedirs(OUT, exist_ok=True)

SIG = {
 'DAM_up': ['APOE','SPP1','LGALS3','CST7','LPL','CTSB','CTSD','CTSL','TREM2','AXL','ITGAX','CD9','CLEC7A','GPNMB','FABP5','LILRB4A','CSF1','B2M','FTH1','LGALS1'],
 'IFN_up': ['IFIT3','ISG15','USP18','IRF7','IFIT1','IFI44','OASL2','RSAD2','STAT1','MX1','IFIT2','BST2','IFIT2','OAS1','IFI27','IFI6','IFITM3','XAF1','EIF2AK2','GBP2'],
 'Homeo_down': ['P2RY12','TMEM119','SALL1','HEXB','SELPLG','CX3CR1','SIGLECH','MAFB','CSF1R','TGFBR1'],
 'Chemokine_up': ['CCL2','CCL3','CCL4','CCL5','CCL12','CXCL10','CXCL2','CXCL16','IL6','TNF'],
}
LIBS = ['LINCS_L1000_Chem_Pert_down', 'LINCS_L1000_Chem_Pert_up']

def enrich(genes, lib):
    r = requests.post('https://maayanlab.cloud/Enrichr/addList',
                      files={'list': (None, '\n'.join(genes)),
                             'description': (None, 'oht')}, timeout=60)
    uid = r.json()['userListId']
    time.sleep(1)
    r = requests.get(f'https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={lib}', timeout=60)
    return pd.DataFrame(r.json()[lib], columns=['rank','term','pval','zscore','combined','genes','adj_pval','old_p','old_adj'])

allres = {}
for sig, genes in SIG.items():
    for lib in LIBS:
        df = enrich(genes, lib)
        df['signature'] = sig; df['library'] = lib
        df['compound'] = df.term.str.split('-').str[0].str.strip()
        allres[f'{sig}|{lib}'] = df
        print(f'{sig} vs {lib}: top {df.iloc[0].term} p={df.iloc[0].pval:.2e} adj={df.iloc[0].adj_pval:.2e}', flush=True)
        time.sleep(2)

res = pd.concat(allres, ignore_index=True)
res.to_csv(f'{OUT}/enrichr_lincs_all.csv', index=False)

# 逆转打分: 上调签名应在 DOWN 库富集; 下调签名(homeo)应在 UP 库富集
rows = []
for comp, gdf in res.groupby('compound'):
    score = 0.0; detail = {}
    for sig, want_lib, w in [('DAM_up','down',1.0), ('IFN_up','down',0.8), ('Chemokine_up','down',0.8), ('Homeo_down','up',1.0)]:
        sub = gdf[(gdf.signature==sig) & (gdf.library.str.endswith(want_lib))]
        if len(sub):
            best = sub.iloc[0]
            s = -np.log10(max(best.adj_pval, 1e-300)) * w
            detail[sig] = f'{best.adj_pval:.1e}'
            score += s
    rows.append({'compound': comp, 'reversal_score': score, **detail})
rev = pd.DataFrame(rows).sort_values('reversal_score', ascending=False)
rev.to_csv(f'{OUT}/drug_reversal_ranking.csv', index=False)
print('\n=== TOP 15 逆转候选 ===', flush=True)
print(rev.head(15).to_string(index=False), flush=True)
print('\nDONE', flush=True)
