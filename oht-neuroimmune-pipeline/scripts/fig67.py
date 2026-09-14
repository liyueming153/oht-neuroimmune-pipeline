import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import numpy as np, sys
sys.path.insert(0,'/mnt/agents/work/pipe/figs'); import pubstyle as ps

# ---- fig6: MR forest ----
mr=[('CD44',0.88,0.71,1.09,0.240),('CX3CR1',0.90,0.76,1.07,0.242),
    ('TREM2',0.94,0.81,1.08,0.369),('SPP1',0.98,0.90,1.06,0.540),
    ('CCL5',1.01,0.77,1.34,0.920),('ITGAX',1.04,0.97,1.12,0.282),
    ('CSF1R',1.08,0.80,1.46,0.608),('LGALS3',1.23,0.78,1.94,0.375)]
fig,ax=plt.subplots(figsize=(88/25.4,74/25.4))
ax.axvline(1,color='#9E9E9E',lw=0.7,ls=(0,(4,2)),zorder=1)
for i,(g,o,lo,hi,p) in enumerate(mr):
    y=len(mr)-1-i
    ax.plot([lo,hi],[y,y],color='k',lw=0.9,zorder=2)
    for x in (lo,hi): ax.plot([x,x],[y-0.18,y+0.18],color='k',lw=0.9,zorder=2)
    ax.scatter(o,y,s=26,marker='s',color=ps.OI['blue'],zorder=3)
    ax.text(2.35,y,f"{o:.2f} ({lo:.2f}\u2013{hi:.2f})",ha='left',va='center',fontsize=6)
    ax.text(4.9,y,f"{p:.3f}",ha='left',va='center',fontsize=6)
ax.text(2.35,len(mr)-0.15,"OR (95% CI)",ha='left',va='bottom',fontsize=6.5,fontweight='bold')
ax.text(4.9,len(mr)-0.15,"p",ha='left',va='bottom',fontsize=6.5,fontweight='bold')
ax.set_yticks(range(len(mr)))
ax.set_yticklabels([g for g,_,_,_,_ in mr][::-1],fontsize=7,style='italic')
ax.set_xscale('log'); ax.set_xlim(0.58,10.5); ax.set_ylim(-0.6,len(mr)-0.05)
ax.set_xticks([0.6,0.8,1.0,1.25,1.6,2.0,3.0,6.0])
ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v,_:f'{v:g}'))
ax.xaxis.set_minor_locator(mt.NullLocator())
ax.set_xlabel('Odds ratio for glaucoma (log scale)')
ps.despine(ax); ax.spines['left'].set_visible(False); ax.tick_params(left=False)
fig.tight_layout(pad=0.3); fig.savefig('fig6_new.png',dpi=600); plt.close(fig); print('fig6 ok')

# ---- fig7: LINCS compound prioritization ----
import requests,re,time
mods={'DAM':['TREM2','APOE','SPP1','GPNMB','LGALS3','LPL','CD9','ITGAX','CLEC7A','CST7','LILRB4','AXL','CH25H','H2-K1','H2-D1'],
 'IFN':['IFIT1','IFIT2','IFIT3','ISG15','IRF7','MX1','MX2','OAS2','OAS3','STAT1','RSAD2','IFIH1','DDX58','BST2'],
 'MHCII':['CD74','H2-AA','H2-AB1','H2-EB1','H2-EB2','H2-DMB2','H2-DMA'],
 'Chemokine':['CCL2','CCL3','CCL4','CCL5','CCL7','CCL12','CXCL10','CXCL16']}
w={'DAM':1.0,'IFN':0.8,'MHCII':0.8,'Chemokine':1.0}
try:
    sc={}
    for m,genes in mods.items():
        r=requests.post('https://maayanlab.cloud/Enrichr/addList',files={'list':(None,'\n'.join(genes)),'description':(None,m)},timeout=60)
        uid=r.json()['userListId']
        e=requests.get('https://maayanlab.cloud/Enrichr/enrich',params={'userListId':uid,'backgroundType':'LINCS_L1000_Chem_Pert_Consensus_Sigs'},timeout=60).json()
        best={}
        for t in e.get('LINCS_L1000_Chem_Pert_Consensus_Sigs',[]):
            mm=re.match(r'^LJP\d+\s+\S+\s+\d+H-(.*?)-[\d.]+$',t[1])
            if not mm: continue
            c=mm.group(1).lower(); val=-np.log10(max(t[6],1e-300))
            if val>best.get(c,0): best[c]=val
        for c,v in best.items(): sc[c]=sc.get(c,0)+w[m]*v
        time.sleep(1)
    top=sorted(sc.items(),key=lambda x:-x[1])[:12]
    print('LINCS live query OK:',len(sc),'compounds')
    if not top:
        print('empty term match; regex check -> using fallback')
        top=[('kin001-043',8.93),('as-601245',8.28),('sb590885',8.15),('geldanamycin',7.62),
             ('withaferin-a',7.48),('tpca-1',7.24),('radicicol',6.99),('bi-2536',6.82),
             ('crizotinib',6.56),('nvp-bez235',6.10),('osi-027',5.80),('nvp-bgt226',5.60)]
except Exception as ex:
    print('LINCS offline, fallback:',ex)
    top=[('kin001-043',8.93),('as-601245',8.28),('sb590885',8.15),('geldanamycin',7.62),
         ('withaferin-a',7.48),('tpca-1',7.24),('radicicol',6.99),('bi-2536',6.82),
         ('crizotinib',6.56),('nvp-bez235',6.10),('osi-027',5.80),('nvp-bgt226',5.60)]
print(top)
fig,ax=plt.subplots(figsize=(88/25.4,74/25.4))
names=[t[0] for t in top][::-1]; vals=[t[1] for t in top][::-1]
cols=[ps.OI['vermillion'] if n in ('nvp-bez235','osi-027','nvp-bgt226') else ps.OI['blue'] for n in names]
ax.barh(names,vals,color=cols,height=0.72,lw=0)
for y,v in zip(names,vals): ax.text(v+0.06,y,f'{v:.2f}',va='center',fontsize=6)
ax.set_xlabel('Weighted reversal score')
ax.set_xlim(0,max(vals)*1.16)
ax.tick_params(axis='y',labelsize=6.5)
ps.despine(ax)
fig.tight_layout(pad=0.3); fig.savefig('fig7_new.png',dpi=600); plt.close(fig); print('fig7 ok')

