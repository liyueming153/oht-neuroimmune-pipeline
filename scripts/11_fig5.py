"""阶段5 Part3: Fig5 绘图 (纯CSV, 无需h5ad)
a: Spp1|Cd44 互作热图 (source x target, Sham vs MB Retina)
b: 关键基因表达点图 (细胞类型 x 基因)
c: GSE304931 签名FC: WT-IOP vs KO-IOP (Spp1依赖性)
d: GSE304931 关键基因表达组间条形图 (Cx3cl1/Cx3cr1/Spp1)
"""
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/tmp/results/05_cci'
plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.8, 'figure.dpi': 150})

# ---------- Fig5a: Spp1|Cd44 互作热图 ----------
focus = pd.read_csv(f'{OUT}/liana_focus_axes.csv')
sp = focus[(focus.axis == 'Spp1|Cd44') & (focus.tissue == 'Retina')].copy()
cts = ['Homeostatic MG','MHCII+ MG','Activated MG (Ccl5+)','DAM-like MG','Iron+ MG (Fth1+)',
       'Proliferating MG','BAM','MHCII+ Mac/APC','pDC','T cell','gdT cell','NK cell','B cell','Neutrophil']
cts = [c for c in cts if c in set(sp.source) | set(sp.target)]
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
for ax, cond in zip(axes, ['Sham', 'MB']):
    d = sp[sp.condition == cond]
    M = pd.DataFrame(0.0, index=cts, columns=cts)
    for _, r in d.iterrows():
        if r.source in cts and r.target in cts:
            M.loc[r.source, r.target] = r.lr_means
    im = ax.imshow(M.values, cmap='Reds', vmin=0, vmax=1.1, aspect='auto')
    ax.set_xticks(range(len(cts)), cts, rotation=90, fontsize=7)
    ax.set_yticks(range(len(cts)), cts, fontsize=7)
    ax.set_title(f'Spp1→CD44 ({cond}, Retina)\nsender=row, receiver=col', fontsize=9)
    ax.set_xlabel('receiver'); ax.set_ylabel('sender' if cond=='Sham' else '')
    # 标注显著格 (magnitude_rank < 0.1)
    for _, r in d.iterrows():
        if r.source in cts and r.target in cts and r.magnitude_rank < 0.1:
            ax.plot(cts.index(r.target), cts.index(r.source), marker='*', color='black', markersize=7)
fig.subplots_adjust(right=0.86)
cb = fig.add_axes([0.88, 0.15, 0.015, 0.7]); fig.colorbar(im, cax=cb, label='LR means')
plt.tight_layout(); plt.savefig(f'{OUT}/Fig5a_Spp1_Cd44_heatmap.png', dpi=300, bbox_inches='tight'); plt.close()
print('Fig5a done', flush=True)

# ---------- Fig5b: 关键基因表达点图 ----------
expr = pd.read_csv(f'{OUT}/key_gene_expr_by_celltype.csv')
genes = ['Cx3cr1','P2ry12','Spp1','Cd44','Ccr5','Cxcr3','C3ar1','Trem2','Axl','Itgav','Itgb1']
genes = [g for g in genes if g in expr.columns]
ct_order = [c for c in cts if c in expr.celltype.values]
fig, ax = plt.subplots(figsize=(7.5, 4.8))
xpos = np.arange(len(genes))
for i, ct in enumerate(ct_order):
    for cond, marker, edge in [('Sham','o','gray'), ('MB','s','red')]:
        sub = expr[(expr.celltype==ct) & (expr.condition==cond)]
        if sub.empty: continue
        for j, g in enumerate(genes):
            pct = sub[g+'_pct'].values[0]; val = sub[g].values[0]
            ax.scatter(j + (0.18 if cond=='MB' else -0.18), i, s=max(pct,0.01)*260,
                       c=val, cmap='viridis', vmin=0, vmax=3.5, marker=marker,
                       edgecolors=edge, linewidths=0.6, alpha=0.9)
ax.set_xticks(xpos, genes, rotation=45, ha='right')
ax.set_yticks(range(len(ct_order)), ct_order, fontsize=8)
ax.set_xlim(-0.6, len(genes)-0.4)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0], marker='o', color='w', markeredgecolor='gray', markersize=8, label='Sham'),
                   Line2D([0],[0], marker='s', color='w', markeredgecolor='red', markersize=8, label='MB (OHT)')],
          loc='upper right', fontsize=8)
ax.set_title('Key ligand/receptor expression by cell type (Retina)\nsize=%expr cells, color=mean log-norm expr; □=MB ○=Sham', fontsize=9)
plt.tight_layout(); plt.savefig(f'{OUT}/Fig5b_keygene_dotplot.png', dpi=300, bbox_inches='tight'); plt.close()
print('Fig5b done', flush=True)

# ---------- Fig5c: GSE304931 签名依赖性 ----------
sd = pd.read_csv(f'{OUT}/GSE304931_signature_dependence.csv')
fig, ax = plt.subplots(figsize=(6.5, 4))
x = np.arange(len(sd)); w = 0.35
ax.bar(x-w/2, sd.WT_IOP_meanFC, w, label='WT + IOP (vs WT)', color='#4C78A8')
ax.bar(x+w/2, sd.KO_IOP_meanFC, w, label='Spp1-KO + IOP (vs KO)', color='#F58518')
for i, r in sd.iterrows():
    for off, v, p in [(-w/2, r.WT_IOP_meanFC, r.WT_IOP_paired_p), (w/2, r.KO_IOP_meanFC, r.KO_IOP_paired_p)]:
        star = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
        ax.text(i+off, v + (0.06 if v>=0 else -0.14), star, ha='center', fontsize=7)
    dp = r.Spp1_dependence_p
    star = '***' if dp<0.001 else '**' if dp<0.01 else '*' if dp<0.05 else 'ns'
    ax.text(i, max(r.WT_IOP_meanFC, r.KO_IOP_meanFC)+0.3, f'dep:{star}', ha='center', fontsize=7, color='crimson')
ax.axhline(0, color='k', lw=0.8)
ax.set_ylim(top=max(sd.WT_IOP_meanFC.max(), sd.KO_IOP_meanFC.max())*1.55, bottom=-0.65)
ax.set_xticks(x, sd.signature, rotation=20, ha='right')
ax.set_ylabel('mean log2FC of signature genes')
ax.set_title('GSE304931 sorted microglia: signature activation with/without Spp1\n(Wilcoxon paired vs 0 above bars; red = WT-vs-KO dependence test)', fontsize=9)
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f'{OUT}/Fig5c_GSE304931_signature_dependence.png', dpi=300, bbox_inches='tight'); plt.close()
print('Fig5c done', flush=True)

# ---------- Fig5d: GSE304931 关键基因分组条形图 ----------
kg = pd.read_csv(f'{OUT}/GSE304931_keygenes.csv')
sel = kg[kg.gene.isin(['Cx3cl1','Cx3cr1','Spp1','Lgals3','Cst7','Ccl5','Cd74','Ifit3'])]
fig, axes = plt.subplots(2, 4, figsize=(11, 5.2))
for ax, (_, r) in zip(axes.ravel(), sel.iterrows()):
    vals = [r.mean_WT, r.mean_IOP, r.mean_KO, r.mean_KO_IOP]
    cols = ['#999999','#4C78A8','#CCCCCC','#F58518']
    ax.bar(['WT','WT+IOP','KO','KO+IOP'], vals, color=cols)
    ax.set_title(r.gene, fontsize=10, style='italic')
    p1 = r.get('IOP_vs_WT_p', np.nan); p2 = r.get('KO_IOP_vs_WT_IOP_p', np.nan)
    lbl = []
    if pd.notna(p1): lbl.append(f'IOPvsWT p={p1:.3f}')
    if pd.notna(p2): lbl.append(f'KOIOPvsWTIOP p={p2:.4f}')
    ax.text(0.5, 0.95, '\n'.join(lbl), transform=ax.transAxes, ha='center', va='top', fontsize=6.5)
    ax.set_ylabel('log2(norm+1)')
fig.suptitle('GSE304931 sorted microglia: brake loss (Cx3cl1/Cx3cr1) & DAM axis genes', fontsize=10)
plt.tight_layout(); plt.savefig(f'{OUT}/Fig5d_GSE304931_keygenes.png', dpi=300, bbox_inches='tight'); plt.close()
print('Fig5d done', flush=True)
print('ALL DONE', flush=True)
