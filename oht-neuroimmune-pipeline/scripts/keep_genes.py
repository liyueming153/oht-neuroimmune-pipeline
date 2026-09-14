# Genes to retain in slim checkpoints beyond the 2500 HVGs:
# annotation signatures (02_annotate), gates, contaminant markers,
# figure markers (fig1c, fig5b), module genes, pseudotime canonicals.
EXTRA_GENES = [
 # ZSIG signatures
 'P2ry12','Tmem119','Cx3cr1','Sall1','Siglech','Hexb','Fcrls','Olfml3','Gpr34','P2ry13',
 'Apoe','Spp1','Lgals3','Cst7','Trem2','Tyrobp','Gpnmb','Itgax','Cd9','Ctsb',
 'Fth1','Ftl1','Ltf','Hmox1','Slc40a1',
 'Ccl5','Ccl3','Ccl4','Cxcl10','Ccl2','Ccl12',
 'Mki67','Top2a','Birc5','Cdk1',
 'Mrc1','Lyve1','Pf4','Folr2','Cd163',
 'Cd74','H2-Aa','H2-Ab1','H2-Eb1','Ciita',
 'Ly6c2','Ccr2','Plac8','Chil3',
 'Bst2','Cox6a2','Klk1','Tcf4',
 'Cd3e','Cd3g','Trbc1','Trbc2',
 'Trdc','Trgc1','Trgc2',
 'Nkg7','Klrb1c','Ncr1','Gzma',
 'Cd79a','Cd79b','Ms4a1','Cd19',
 'S100a8','S100a9','Ly6g','Cxcr2','Mmp9','Mpo',
 # contaminant markers
 'Rho','Pde6a','Pde6b','Nrl','Opn1sw','Rbpms','Nefl','Sncg','Slc17a6',
 'Rlbp1','Glul','Slc1a3','Clu','Gfap','Aqp4','S100b','Aldh1l1',
 'Pecam1','Kdr','Cldn5','Rpe65','Mitf','Tyr',
 # figure markers / fig5b / canonicals
 'Cd44','Ccr5','Cxcr3','C3ar1','Axl','Itgav','Itgb1','Lyz2','Cd63','Ctsd',
 'Ifit1','Ifit2','Ifit3','Isg15','Irf7','Mx1','Mx2','Oas2','Oas3','Stat1','Rsad2','Usp18','Ifih1','Ddx58',
 'H2-Eb2','H2-DMb2','H2-DMa','H2-D1','H2-K1',
 'Ccl7','Cxcl16','Cxcl2','Il1b','Nfkb1',
 'Cx3cl1','Oasl2','Csf1','Rpl3','Rplp0','Rps3','Rps18','Rps27','Flt1','Gzmb','Prf1','Xcl1','Klrd1','Il7r','Il23r','Ltb','Icos',
]
