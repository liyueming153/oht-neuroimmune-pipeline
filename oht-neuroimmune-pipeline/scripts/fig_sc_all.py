#!/usr/bin/env python3
"""Single-cell figures 1a-d, 2a-c, 3a-d, 5b in publication style + stats CSVs."""
import sys, os
sys.path.insert(0, '/mnt/agents/work/pipe/figs')
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
import matplotlib.pyplot as plt
import pubstyle as ps
ps.apply()
OUT = os.environ.get('OUT_DIR', '/tmp/ohtrun/out')

ad = sc.read_h5ad(f'{OUT}/annotated.h5ad')
ad = ad[ad.obs.cell_class != 'Contaminant']
ad.obs['grp'] = ad.obs.tissue.astype(str) + '-' + ad.obs.condition.astype(str)
CLASSES = ['Homeostatic MG','DAM-like MG','Activated MG (Ccl5+)','Proliferating MG',
           'MHCII+ Mac/APC','pDC','B cell','T cell','gdT cell','NK cell','Neutrophil']
present = [c for c in CLASSES if c in set(ad.obs.cell_class)]
xy = ad.obsm['X_umap']

# ---------- Fig 1a: UMAP by class
fig, ax = plt.subplots(figsize=(4.6, 3.2))
ps.umap_panel(ax, xy, ad.obs.cell_class.values, ps.CELL, s=0.5, legend=True)
ax.set_title('GSE293358 CD45+ immune cells')
fig.savefig(f'{OUT}/fig1a_new.png'); plt.close(fig); print('1a ok', flush=True)

# ---------- Fig 1b: UMAP by condition
fig, ax = plt.subplots(figsize=(4.0, 3.2))
cond = ad.obs.condition.values
ps.umap_panel(ax, xy, cond, {'Sham': ps.OI['grey'], 'MB': ps.OI['vermillion']}, s=0.5, legend=True)
ax.set_title('Condition')
fig.savefig(f'{OUT}/fig1b_new.png'); plt.close(fig); print('1b ok', flush=True)

# ---------- Fig 1c: marker dotplot
MARKERS = ['P2ry12','Tmem119','Sall1','Apoe','Spp1','Lgals3','Cst7','Fth1','Ftl1',
           'Mki67','Top2a','Cd74','H2-Aa','H2-Ab1','Ly6c2','Ccr2','Cd79a','Ms4a1',
           'Cd3e','Cd3g','Trdc','Nkg7','Klrb1c','S100a8','S100a9','Ly6g']
MARKERS = [g for g in MARKERS if g in ad.var_names]
ad.obs['cell_class'] = pd.Categorical(ad.obs['cell_class'], categories=present, ordered=True)
dp = sc.pl.DotPlot(ad, MARKERS, groupby='cell_class')
mean = dp.dot_color_df.loc[present, MARKERS]
frac = dp.dot_size_df.loc[present, MARKERS]
fig, ax = plt.subplots(figsize=(6.4, 3.0))
xmax = frac.values.max()
for i, c in enumerate(present):
    for j, g in enumerate(MARKERS):
        ax.scatter(j, i, s=28*frac.loc[c, g]/xmax if xmax else 1,
                   c=mean.loc[c, g], cmap='viridis', vmin=0, vmax=mean.values.max(),
                   linewidths=0, )
ax.set_xticks(range(len(MARKERS))); ax.set_xticklabels(MARKERS, rotation=45, ha='right', style='italic')
ax.set_yticks(range(len(present))); ax.set_yticklabels(present)
ax.set_xlim(-0.6, len(MARKERS)-0.4); ax.set_ylim(len(present)-0.5, -0.6)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, mean.values.max()))
cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.12); cb.set_label('Mean expr.', fontsize=6)
cb.outline.set_visible(False)
# size legend
for k, f in enumerate([0.25, 0.5, 0.75, 1.0]):
    ax.scatter(len(MARKERS)+0.4+k*0.9, len(present)-0.2, s=28*f, c='grey', linewidths=0)
    ax.text(len(MARKERS)+0.4+k*0.9, len(present)+0.55, f'{int(f*100)}', ha='center', fontsize=5.5)
ax.text(len(MARKERS)+1.6, len(present)+1.3, '% expr.', fontsize=6, ha='center')
ax.set_xlim(-0.6, len(MARKERS)+3.6)
fig.savefig(f'{OUT}/fig1c_new.png'); plt.close(fig); print('1c ok', flush=True)

# ---------- Fig 1d: composition stacked bars + per-sample stats
GRPS = ['OpticNerve-Sham','OpticNerve-MB','Retina-Sham','Retina-MB']
tab = pd.crosstab(ad.obs.grp, ad.obs.cell_class)[present].loc[GRPS]
prop = tab.div(tab.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(3.4, 3.0))
bottom = np.zeros(4)
for c in present:
    ax.bar(range(4), prop[c], bottom=bottom, color=ps.CELL.get(c,'#BBBBBB'),
           width=0.62, label=c, linewidth=0)
    bottom += prop[c].values
ax.set_xticks(range(4)); ax.set_xticklabels(['ON\nSham','ON\nMB','Retina\nSham','Retina\nMB'])
ax.set_ylabel('Fraction of CD45+ cells'); ax.set_ylim(0,1)
ax.legend(loc='center left', bbox_to_anchor=(1.02,0.5), handlelength=1.0, labelspacing=0.3)
ps.despine(ax)
fig.savefig(f'{OUT}/fig1d_new.png'); plt.close(fig); print('1d ok', flush=True)

# per-sample composition stats (retina + ON), MB vs Sham, Welch on log2 prop
rows=[]
sam = pd.crosstab([ad.obs.tissue, ad.obs.condition, ad.obs['sample']], ad.obs.cell_class)
sam = sam[sam.sum(axis=1)>0]
sprop = sam.div(sam.sum(axis=1), axis=0)
for tissue in ['Retina','OpticNerve']:
    for c in present:
        a = np.log2(sprop.loc[(tissue,'MB')][c]+1e-4)
        b = np.log2(sprop.loc[(tissue,'Sham')][c]+1e-4)
        t,p = stats.ttest_ind(a,b,equal_var=False)
        rows.append(dict(tissue=tissue, cell_class=c, log2fc=a.mean()-b.mean(), p=p,
                         prop_MB=sprop.loc[(tissue,'MB')][c].mean(), prop_Sham=sprop.loc[(tissue,'Sham')][c].mean()))
st = pd.DataFrame(rows)
for tissue in ['Retina','OpticNerve']:
    m = st.tissue==tissue
    st.loc[m,'fdr'] = stats.false_discovery_control(st.loc[m,'p'].fillna(1))
st.to_csv(f'{OUT}/composition_stats.csv', index=False)
print(st[(st.tissue=='Retina')].sort_values('p').to_string(), flush=True)

# retina microglial homeostatic fraction of MG
mgmask = ad.obs.cell_class.isin(['Homeostatic MG','DAM-like MG','Iron+ MG (Fth1+)','Proliferating MG'])
mg = ad.obs[mgmask]
for cond in ['Sham','MB']:
    m2 = (mg.tissue=='Retina')&(mg.condition==cond)
    frac_h = (mg[m2].cell_class=='Homeostatic MG').mean()
    print(f'Retina {cond}: homeostatic fraction of MG = {frac_h:.3f}, nMG={m2.sum()}', flush=True)
del ad
import gc; gc.collect()
print('PART1 DONE', flush=True)

# ================= PART 2: microglia figures 2a-c, 3a-d =================
mg = sc.read_h5ad(f'{OUT}/microglia_dpt.h5ad')
mg.obs['grp'] = mg.obs.tissue.astype(str) + '-' + mg.obs.condition.astype(str)
MODS = ['homeostatic','DAM','MHCII','interferon','chemokine','iron','proliferating','ribosomal']
Zc = mg.obs.groupby('mg_leiden', observed=True)[[f'mod_{m}' for m in MODS]].mean()
Zc.columns = MODS
# label clusters
def label_cluster(r):
    if r['proliferating'] > 0.4 and r['proliferating'] >= r.drop(['proliferating','ribosomal']).max():
        return 'MG-Prolif'
    if r['ribosomal'] > 0.5 and r.drop(['ribosomal']).max() < 0.3:
        return 'MG-lowRNA'
    act = r[['DAM','MHCII','interferon','chemokine','iron']]
    if r['homeostatic'] > act.max() + 0.1:
        return None  # homeostatic; numbered later
    top = act.idxmax()
    return {'DAM':'MG-DAM','MHCII':'MG-MHCII+','interferon':'MG-IFN',
            'chemokine':'MG-Ccl5+','iron':'MG-DAM-Iron'}[top]
lab = {}
homeo_clusters = [c for c in Zc.index if label_cluster(Zc.loc[c]) is None]
hord = Zc.loc[homeo_clusters,'homeostatic'].sort_values(ascending=False).index
for i,c in enumerate(hord,1): lab[c] = 'MG-Homeo' if len(hord)==1 else f'MG-Homeo-{i}'
for c in Zc.index:
    if c not in lab: lab[c] = label_cluster(Zc.loc[c])
mg.obs['mg_label'] = mg.obs['mg_leiden'].map(lab).astype(str)
STATES = [lab[c] for c in Zc['homeostatic'].sort_values(ascending=False).index]
STATES = sorted(set(mg.obs['mg_label']), key=lambda s:(not s.startswith('MG-Homeo'), s))
STATE_C = dict(zip(STATES, ['#1b9e77','#66c2a5','#a6d854','#01665e','#8da0cb',
                            '#e78ac3','#d95f02','#e6ab02','#a6761d','#7570b3','#666666','#999999']))
print('states:', STATES, flush=True)

# ---------- Fig 2a: MG UMAP by state
fig, ax = plt.subplots(figsize=(4.6,3.2))
ps.umap_panel(ax, mg.obsm['X_umap'], mg.obs['mg_label'].values, STATE_C, s=0.5, legend=True)
ax.set_title(f'Microglia re-clustering (n={mg.n_obs:,})')
fig.savefig(f'{OUT}/fig2a_new.png'); plt.close(fig); print('2a ok', flush=True)

# ---------- Fig 2b: group-split UMAP
fig, axes = plt.subplots(2,2, figsize=(5.4,5.0))
for ax, g in zip(axes.ravel(), ['Retina-Sham','Retina-MB','OpticNerve-Sham','OpticNerve-MB']):
    m = (mg.obs.grp==g).values
    ax.scatter(mg.obsm['X_umap'][:,0], mg.obsm['X_umap'][:,1], s=0.4, c='#E0E0E0', linewidths=0, rasterized=True)
    col = ps.OI['vermillion'] if g.endswith('MB') else '#5E9C94'
    ax.scatter(mg.obsm['X_umap'][m,0], mg.obsm['X_umap'][m,1], s=0.4, c=col, linewidths=0, rasterized=True)
    ax.set_title(f'{g} (n={m.sum():,})')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
fig.tight_layout()
fig.savefig(f'{OUT}/fig2b_new.png'); plt.close(fig); print('2b ok', flush=True)

# ---------- Fig 2c: MG-state composition bars + stats
GRPS = ['OpticNerve-Sham','OpticNerve-MB','Retina-Sham','Retina-MB']
tab = pd.crosstab(mg.obs.grp, mg.obs.mg_label).reindex(GRPS).fillna(0)
tab = tab[[c for c in STATES if c in tab.columns]]
prop = tab.div(tab.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(3.6,3.0))
bottom = np.zeros(4)
for c in tab.columns:
    ax.bar(range(4), prop[c], bottom=bottom, color=STATE_C[c], width=0.62, label=c, linewidth=0)
    bottom += prop[c].values
ax.set_xticks(range(4)); ax.set_xticklabels(['ON\nSham','ON\nMB','Retina\nSham','Retina\nMB'])
ax.set_ylabel('Fraction of microglia'); ax.set_ylim(0,1)
ax.legend(loc='center left', bbox_to_anchor=(1.02,0.5), handlelength=1.0, labelspacing=0.3)
ps.despine(ax)
fig.savefig(f'{OUT}/fig2c_new.png'); plt.close(fig); print('2c ok', flush=True)

rows=[]
sam = pd.crosstab([mg.obs.tissue, mg.obs.condition, mg.obs['sample']], mg.obs.mg_label)
sam = sam[sam.sum(axis=1)>0]
sprop = sam.div(sam.sum(axis=1), axis=0)
for tissue in ['Retina','OpticNerve']:
    for c in tab.columns:
        a = np.log2(sprop.loc[(tissue,'MB')][c]+1e-4) if (tissue,'MB') in sprop.index else None
        b = np.log2(sprop.loc[(tissue,'Sham')][c]+1e-4)
        if a is None or len(a)<2: continue
        t,p = stats.ttest_ind(a,b,equal_var=False)
        rows.append(dict(tissue=tissue, state=c, log2fc=a.mean()-b.mean(), p=p))
mst = pd.DataFrame(rows)
mst.to_csv(f'{OUT}/mgstate_stats.csv', index=False)
print(mst[mst.tissue=='Retina'].sort_values('log2fc').to_string(), flush=True)

# ================= PART 3: pseudotime figures 3a-d =================
# ---------- Fig 3a: DPT UMAP
fig, ax = plt.subplots(figsize=(3.6,3.0))
dpt = mg.obs['dpt_pseudotime'].values
sc_ = ax.scatter(mg.obsm['X_umap'][:,0], mg.obsm['X_umap'][:,1], s=0.5, c=dpt,
                 cmap='magma', linewidths=0, rasterized=True)
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(sc_, ax=ax, fraction=0.04, pad=0.02)
cb.set_label('DPT pseudotime', fontsize=6); cb.outline.set_visible(False)
ax.set_title('Diffusion pseudotime')
fig.savefig(f'{OUT}/fig3a_new.png'); plt.close(fig); print('3a ok', flush=True)

# ---------- Fig 3b: density by group
from scipy.stats import gaussian_kde, mannwhitneyu
fig, ax = plt.subplots(figsize=(3.2,2.6))
GC = {'Retina-Sham':'#5E9C94','Retina-MB':ps.OI['vermillion'],
      'OpticNerve-Sham':'#8C8C8C','OpticNerve-MB':ps.OI['orange']}
xs = np.linspace(0,1,300)
for g,c in GC.items():
    v = dpt[(mg.obs.grp==g).values]
    v = v[np.isfinite(v)]
    kde = gaussian_kde(v, bw_method=0.06)
    ax.plot(xs, kde(xs), color=c, lw=1.2, label=f'{g} (n={len(v):,})')
ax.set_xlabel('Pseudotime (DPT)'); ax.set_ylabel('Density')
ax.legend(handlelength=1.2, labelspacing=0.3)
ps.despine(ax)
fig.tight_layout()
fig.savefig(f'{OUT}/fig3b_new.png'); plt.close(fig); print('3b ok', flush=True)

# ---------- Fig 3c: boxplots per tissue + MWU
fig, axes = plt.subplots(1,2, figsize=(4.0,2.4), sharey=True)
stat_rows=[]
for ax, tissue in zip(axes, ['Retina','OpticNerve']):
    a = dpt[((mg.obs.tissue==tissue)&(mg.obs.condition=='Sham')).values]
    b = dpt[((mg.obs.tissue==tissue)&(mg.obs.condition=='MB')).values]
    a,b = a[np.isfinite(a)], b[np.isfinite(b)]
    bp = ax.boxplot([a,b], widths=0.5, showfliers=False, patch_artist=True,
                    medianprops=dict(color='black', lw=0.8),
                    boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    for patch, c in zip(bp['boxes'], ['#8C8C8C', ps.OI['vermillion']]):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    u,p = mannwhitneyu(a,b)
    dm = np.median(b)-np.median(a)
    stat_rows.append(dict(tissue=tissue, dmedian=dm, mwu_p=p, n_sham=len(a), n_mb=len(b)))
    ax.set_title(f'{tissue} (Δmedian={dm:+.4f})', fontsize=6.5)
    ax.set_xticklabels(['Sham','MB'])
    ps.despine(ax)
axes[0].set_ylabel('Pseudotime')
fig.text(0.5,0.005,'cell-level comparison, descriptive', ha='center', fontsize=6)
fig.tight_layout()
fig.savefig(f'{OUT}/fig3c_new.png'); plt.close(fig)
pd.DataFrame(stat_rows).to_csv(f'{OUT}/dpt_group_tests.csv', index=False)
print('3c ok', stat_rows, flush=True)

# ---------- Fig 3d: pseudotime-binned heatmap
traj = pd.read_csv(f'{OUT}/dpt_trajectory_genes.csv')
canon = ['P2ry12','Tmem119','Sall1','Cx3cr1','Apoe','Ctsb','Ctsd','Cd63','Lyz2','Trem2',
         'Spp1','Cst7','Cd74','H2-Aa','Ifit3','Isg15','Ccl5','Lgals3','Fth1','Mki67','Usp18']
top_pos = traj.sort_values('rho',ascending=False).head(10).gene.tolist()
top_neg = traj.sort_values('rho').head(10).gene.tolist()
genes = list(dict.fromkeys(top_neg + canon + top_pos))
genes = [g for g in genes if g in mg.var_names][:45]
rng = np.random.default_rng(7)
take = rng.choice(mg.n_obs, size=min(8000, mg.n_obs), replace=False)
sub = mg[take]
d = sub.obs['dpt_pseudotime'].values
ok = np.isfinite(d)
sub = sub[ok]; d = d[ok]
bins = np.linspace(0,1,51)
ib = np.digitize(d, bins)-1
X = sub[:, genes].X
X = X.toarray() if hasattr(X,'toarray') else X
M = np.full((len(genes), 50), np.nan)
for k in range(50):
    m = ib==k
    if m.sum()>0: M[:,k] = X[m].mean(axis=0)
Z = (M - np.nanmean(M,axis=1,keepdims=True)) / (np.nanstd(M,axis=1,keepdims=True)+1e-9)
fig, ax = plt.subplots(figsize=(3.4,4.6))
im = ax.imshow(Z, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, style='italic', fontsize=5)
ax.set_xticks([0,12,25,37,49]); ax.set_xticklabels(['0.0','0.25','0.5','0.75','1.0'])
ax.set_xlabel('Pseudotime (homeostatic → activated)')
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label('z-scored expr.', fontsize=6); cb.outline.set_visible(False)
fig.savefig(f'{OUT}/fig3d_new.png'); plt.close(fig); print('3d ok', flush=True)

# ================= PART 4: Fig 5b L/R dotplot (retina) =================
ad = sc.read_h5ad(f'{OUT}/annotated.h5ad')
ad = ad[(ad.obs.cell_class!='Contaminant') & (ad.obs.tissue=='Retina')]
CLASSES = ['Homeostatic MG','DAM-like MG','Activated MG (Ccl5+)','Proliferating MG',
           'MHCII+ Mac/APC','pDC','B cell','T cell','gdT cell','NK cell','Neutrophil']
present = [c for c in CLASSES if c in set(ad.obs.cell_class)]
GENES = ['Cx3cr1','Spp1','Cd44','Ccr5','Cxcr3','C3ar1','Trem2','Axl','Itgav','Itgb1']
GENES = [g for g in GENES if g in ad.var_names]
X = ad[:, GENES].X
X = X.toarray() if hasattr(X,'toarray') else np.asarray(X)
dfx = pd.DataFrame(X, columns=GENES)
dfx['cls'] = ad.obs.cell_class.values
dfx['cond'] = ad.obs.condition.values
fig, ax = plt.subplots(figsize=(5.6,3.2))
vmax = 0
cells = {}
for i,cname in enumerate(present):
    for j,g in enumerate(GENES):
        for k,(cond,col) in enumerate([('Sham','#8C8C8C'),('MB',ps.OI['vermillion'])]):
            v = dfx[(dfx.cls==cname)&(dfx.cond==cond)][g]
            if len(v)==0: continue
            fr = (v>0).mean(); me = v[v>0].mean() if (v>0).any() else 0
            vmax = max(vmax, me)
            cells[(i,j,k)] = (fr, me)
norm = plt.Normalize(0, vmax)
for (i,j,k),(fr,me) in cells.items():
    ax.scatter(j + (k-0.5)*0.34, i, s=90*fr, c=[plt.cm.viridis(norm(me))],
               linewidths=0.5, edgecolors='#8C8C8C' if k==0 else ps.OI['vermillion'])
ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=45, ha='right', style='italic')
ax.set_yticks(range(len(present))); ax.set_yticklabels(present)
ax.set_xlim(-0.6, len(GENES)-0.4); ax.set_ylim(len(present)-0.6, -0.6)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.10); cb.set_label('Mean expr. (expr. cells)', fontsize=6)
cb.outline.set_visible(False)
hs = [plt.Line2D([],[], marker='o', ls='', mfc='white', mec=c, ms=6, label=l)
      for l,c in [('Sham','#8C8C8C'),('MB (OHT)',ps.OI['vermillion'])]]
ax.legend(handles=hs, loc='lower left', bbox_to_anchor=(1.03, 0.0), handletextpad=0.2)
fig.savefig(f'{OUT}/fig5b_new.png'); plt.close(fig); print('5b ok', flush=True)
print('PARTS 2-4 DONE', flush=True)
