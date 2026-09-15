"""导出全部图版的 GraphPad Prism 作图原始数据 → results/plot_data/"""
import os, gc, subprocess
import numpy as np, pandas as pd
from scipy import stats

BASE = "/mnt/agents/output/oht_project"
OUTD = f"{BASE}/results/plot_data"
os.makedirs(OUTD, exist_ok=True)

def log(*a):
    print(*a, flush=True)

def rss():
    with open('/proc/self/status') as f:
        for l in f:
            if l.startswith('VmRSS'): return float(l.split()[1])/1e6

# ============ 1) GSE293358_final.h5ad：Fig1a/1b/1c/1d ============
import scanpy as sc
import scipy.sparse as sp
if os.path.exists(f"{OUTD}/Fig1d_细胞组成百分比.csv"):
    ad = None
    log("Fig1 系列已存在，跳过")
else:
    subprocess.run(f"cat {BASE}/data/GSE293358_final.h5ad.part_* > /tmp/final.h5ad", shell=True, check=True)
    ad = sc.read_h5ad("/tmp/final.h5ad")
if ad is not None:
    log(f"final.h5ad {ad.shape} RSS {rss():.2f}")

    # Fig1a/1b: UMAP 坐标 + 注释
    CT_ORDER = ["Homeostatic MG","MHCII+ MG","Activated MG (Ccl5+)","DAM-like MG",
                "Iron+ MG (Fth1+)","Proliferating MG","BAM","MHCII+ Mac/APC",
                "pDC","T cell","gdT cell","NK cell","B cell","Neutrophil"]
    um = pd.DataFrame(ad.obsm["X_umap"], columns=["UMAP1","UMAP2"], index=ad.obs_names)
    meta = ad.obs[["celltype","condition","tissue","sample"]].copy()
    um = pd.concat([um, meta], axis=1)
    um.to_csv(f"{OUTD}/Fig1ab_UMAP坐标_全量.csv")
    # 分层抽样 ~12000 供 Prism
    ds = um.groupby("celltype", group_keys=False).apply(
        lambda g: g.sample(frac=min(1.0, 12000/len(um)), random_state=0))
    ds.to_csv(f"{OUTD}/Fig1ab_UMAP坐标_精简版.csv", index=False)
    log(f"UMAP 全量 {len(um)} 精简 {len(ds)}")

    # Fig1c: dotplot 统计（mean expr + pct expressing）
    MK = {"Homeostatic MG":["P2ry12","Tmem119","Sall1"],"MHCII+ MG":["H2-Aa","H2-Ab1"],
          "Activated MG (Ccl5+)":["Ccl5","Apoe"],"DAM-like MG":["Cd63","Ctsb","Ctsd"],
          "Iron+ MG (Fth1+)":["Fth1","Ftl1"],"Proliferating MG":["Mki67","Top2a"],
          "BAM":["Mrc1","Ms4a7"],"MHCII+ Mac/APC":["Cd74","H2-Eb1"],"pDC":["Tcf4","Bst2"],
          "T cell":["Cd3e","Cd8a"],"gdT cell":["Trdc","Il7r"],"NK cell":["Nkg7","Ncr1"],
          "B cell":["Cd79a","Ms4a1"],"Neutrophil":["S100a8","S100a9"]}
    mk_flat = [g for v in MK.values() for g in v if g in ad.var_names]
    X = ad[:, mk_flat].X
    if sp.issparse(X): X = X.toarray()
    X = np.asarray(X)
    ct = ad.obs["celltype"].astype(str).to_numpy()
    rows_mean, rows_pct = {}, {}
    for c in CT_ORDER:
        m = ct == c
        sub = X[m]
        rows_mean[c] = sub.mean(axis=0)
        rows_pct[c] = (sub > 0).mean(axis=0) * 100
    pd.DataFrame(rows_mean, index=mk_flat).round(4).to_csv(f"{OUTD}/Fig1c_dotplot_平均表达.csv")
    pd.DataFrame(rows_pct, index=mk_flat).round(2).to_csv(f"{OUTD}/Fig1c_dotplot_表达细胞百分比.csv")
    log("Fig1c done")

    # Fig1d: 组成（4 列 × 14 细胞类型，%）
    grp = pd.crosstab([ad.obs["tissue"], ad.obs["condition"]], ad.obs["celltype"])
    grp_prop = (grp.div(grp.sum(axis=1), axis=0) * 100).T.loc[CT_ORDER]
    grp_prop.to_csv(f"{OUTD}/Fig1d_细胞组成百分比.csv")
    del ad, X, um, ds; gc.collect(); os.remove("/tmp/final.h5ad")
    log(f"final 释放 RSS {rss():.2f}")

# ============ 2) GSE293358_mg.h5ad：Fig2a/2b/2c、Fig3a-d ============
subprocess.run(f"cat {BASE}/data/GSE293358_mg.h5ad.part_* > /tmp/mg.h5ad", shell=True, check=True)
mg = sc.read_h5ad("/tmp/mg.h5ad")
log(f"mg.h5ad {mg.shape} RSS {rss():.2f}")

STATE_LABEL = {"0":"MG-Homeo-ON","1":"MG-Homeo-Ret","2":"MG-Homeo-Ret2","3":"MG-Homeo-ON2",
               "4":"MG-MHCII+","5":"MG-lowRNA","6":"MG-IFN","7":"MG-Ccl5+ act.",
               "8":"MG-Prolif","9":"MG-DAM","10":"MG-MHCII+ (Ret)","11":"MG-DAM-Iron"}
mg.obs["mg_label"] = mg.obs["mg_state"].astype(str).map(STATE_LABEL)
mg.obs["group"] = (mg.obs["tissue"].astype(str) + "-" + mg.obs["condition"].astype(str))
um2 = pd.DataFrame(mg.obsm["X_umap"], columns=["UMAP1","UMAP2"], index=mg.obs_names)
keep = [c for c in ["mg_state","mg_label","group","tissue","condition","sample","dpt_pseudotime"] if c in mg.obs.columns]
um2 = pd.concat([um2, mg.obs[keep]], axis=1)
um2.to_csv(f"{OUTD}/Fig2ab3ac_MG_UMAP坐标与DPT_全量.csv")
ds2 = um2.groupby("mg_label", group_keys=False).apply(
    lambda g: g.sample(frac=min(1.0, 10000/len(um2)), random_state=0))
ds2.to_csv(f"{OUTD}/Fig2ab3ac_MG_UMAP坐标与DPT_精简版.csv", index=False)
log(f"MG UMAP 全量 {len(um2)} 精简 {len(ds2)}")

# Fig3b: KDE 曲线（4 组，与原图同参数）
xs = np.linspace(0, np.nanmax(mg.obs["dpt_pseudotime"]), 300)
kde_df = pd.DataFrame({"pseudotime": xs})
for gname in ["Retina-Sham","Retina-MB","OpticNerve-Sham","OpticNerve-MB"]:
    v = mg.obs.loc[mg.obs["group"] == gname, "dpt_pseudotime"].to_numpy()
    v = v[np.isfinite(v)]
    kde_df[gname] = stats.gaussian_kde(v)(xs)
kde_df.round(5).to_csv(f"{OUTD}/Fig3b_DPT密度曲线.csv", index=False)
# Fig3c: 原始 DPT 值（4 组，分列）
cols = {}
for gname in ["Retina-Sham","Retina-MB","OpticNerve-Sham","OpticNerve-MB"]:
    v = mg.obs.loc[mg.obs["group"] == gname, "dpt_pseudotime"].dropna().to_numpy()
    cols[gname] = pd.Series(v)
pd.DataFrame(cols).round(4).to_csv(f"{OUTD}/Fig3c_DPT原始值.csv", index=False)
log("Fig3b/c done")

# Fig3d: 轨迹热图矩阵（复现 06 脚本：top12负 + 20已知 + top12正，120分箱）
traj = pd.read_csv(f"{BASE}/results/03_mg/traj_gene_correlations.csv")
KNOWN = ["P2ry12","Tmem119","Sall1","Cx3cr1","Apoe","Ctsb","Ctsd","Cd63",
         "Lyz2","Trem2","Spp1","Cst7","Cd74","H2-Aa","Ifit3","Isg15",
         "Ccl5","Lgals3","Fth1","Mki67"]
top_neg = traj.head(12)["gene"].tolist()
top_pos = traj.tail(12)["gene"].tolist()[::-1]
genes_hm = [g for g in dict.fromkeys(top_neg + KNOWN + top_pos) if g in mg.var_names][:48]
dpt = mg.obs["dpt_pseudotime"].to_numpy()
valid = np.isfinite(dpt)
order = np.argsort(dpt[valid])
Xv = mg[valid][:, genes_hm].X
if sp.issparse(Xv): Xv = Xv.toarray()
Xv = np.asarray(Xv)[order]
bins = np.array_split(np.arange(Xv.shape[0]), 120)
Hm = np.vstack([Xv[b].mean(axis=0) for b in bins])
Hm = np.log1p(Hm * 100)
Hm = (Hm - Hm.mean(axis=0)) / (Hm.std(axis=0) + 1e-9)
hm_df = pd.DataFrame(Hm.T, index=genes_hm, columns=[f"bin{i+1}" for i in range(120)])
hm_df.round(3).to_csv(f"{OUTD}/Fig3d_轨迹热图矩阵_zscore.csv")
# 每个分箱的伪时间中点，供横轴标注
bmid = [float(np.sort(dpt[valid])[b].mean()) for b in bins]
pd.DataFrame({"bin": [f"bin{i+1}" for i in range(120)], "pseudotime_mean": bmid}).to_csv(
    f"{OUTD}/Fig3d_分箱伪时间表.csv", index=False)
log(f"Fig3d done, genes={len(genes_hm)}")
del mg, Xv, Hm; gc.collect(); os.remove("/tmp/mg.h5ad")
log(f"mg 释放 RSS {rss():.2f}")

# ============ 3) GSE304931：Fig5d 逐样本表达 ============
d = pd.read_csv(f"{BASE}/data/GSE304931.csv.gz")
d.columns = [str(c) for c in d.columns]
idcol = d.columns[0]
d["ensembl"] = d[idcol].astype(str).str.split(".").str[0]
import mygene
mgq = mygene.MyGeneInfo()
res = mgq.querymany(["Cx3cl1","Cx3cr1","Spp1","Lgals3","Cst7","Ccl5","Cd74","Ifit3"],
                    scopes="symbol", fields="ensembl.gene", species="mouse", as_dataframe=False)
sym2ens = {}
for r in res:
    sym = r.get("query")
    ens = r.get("ensembl")
    if isinstance(ens, dict): ens = [ens]
    if ens: sym2ens[sym] = [e["gene"] for e in ens if isinstance(e, dict) and "gene" in e]
log(sym2ens)
GROUPS = {"WT": ["WT_1","WT_2","WT_3","WT_4"], "WT_IOP": ["IOP_1","IOP_2","IOP_3","IOP_4"],
          "KO": ["Spp1KO_1","Spp1KO_2","Spp1KO_3","Spp1KO_4"],
          "KO_IOP": ["Spp1KO_IOP_1","Spp1KO_IOP_2","Spp1KO_IOP_3","Spp1KO_IOP_4"]}
samp_cols = sum(GROUPS.values(), [])
have = [c for c in samp_cols if c in d.columns]
log(f"样本列命中 {len(have)}/{len(samp_cols)}: {have[:6]}...")
rows = []
for sym, ens_list in sym2ens.items():
    sub = d[d["ensembl"].isin(ens_list)]
    if len(sub) == 0:
        log(f"!! {sym} 未命中"); continue
    vals = np.log2(sub[have].astype(float).mean(axis=0) + 1)  # 同 symbol 多 Ensembl 取均值
    for s in have:
        gname = [k for k, v in GROUPS.items() if s in v][0]
        rows.append({"gene": sym, "sample": s, "group": gname, "log2_expr": round(float(vals[s]), 4)})
pd.DataFrame(rows).to_csv(f"{OUTD}/Fig5d_关键基因逐样本表达.csv", index=False)
log("Fig5d done")
log("ALL EXPORT DONE")
