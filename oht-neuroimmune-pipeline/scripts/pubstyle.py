"""Publication figure style — Nature/Cell-family bioinformatics conventions.
Arial/Helvetica 7pt base, Okabe-Ito colorblind-safe palette, despined 0.7pt
axes, no gridlines, 300 dpi."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

OI = dict(black='#000000', orange='#E69F00', sky='#56B4E9', green='#009E73',
          yellow='#F0E442', blue='#0072B2', vermillion='#D55E00', pink='#CC79A7',
          grey='#8C8C8C', lightgrey='#D9D9D9')
COND = {'Sham': OI['grey'], 'MB': OI['vermillion'], 'OHT': OI['vermillion'],
        'naive': OI['grey'], 'WT': OI['grey'], 'KO': OI['sky']}
CELL = {
 'Homeostatic MG': '#56B4E9', 'DAM-like MG': '#D55E00', 'MHCII+ Mac/APC': '#009E73',
 'Iron+ MG (Fth1+)': '#CC79A7', 'Proliferating MG': '#F0E442', 'BAM': '#0072B2',
 'Activated MG (Ccl5+)': '#E69F00', 'Monocyte': '#8C8C8C', 'Neutrophil': '#000000',
 'B cell': '#88CCEE', 'T cell': '#44AA99', 'NK cell': '#117733', 'gdT cell': '#999933',
 'pDC': '#CC79A7', 'Contaminant': '#D9D9D9', 'Other': '#DDDDDD',
}

def apply():
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
        'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
        'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
        'xtick.major.size': 3, 'ytick.major.size': 3,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False,
        'legend.frameon': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    })

def despine(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3, width=0.7)

def umap_panel(ax, xy, labels, colors, s=0.7, legend=False, title=None, rasterized=True):
    import numpy as np
    labs = np.asarray(labels)
    for lab in dict.fromkeys(labs):
        m = labs == lab
        ax.scatter(xy[m, 0], xy[m, 1], s=s, c=colors.get(lab, '#BBBBBB'),
                   linewidths=0, label=lab, rasterized=rasterized)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    if title: ax.set_title(title)
    if legend:
        ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), markerscale=3,
                  handletextpad=0.1, borderaxespad=0, labelspacing=0.35)
    return ax

def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'
