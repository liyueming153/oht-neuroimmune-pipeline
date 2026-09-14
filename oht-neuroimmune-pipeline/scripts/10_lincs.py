#!/usr/bin/env python3
"""Step 10: LINCS L1000 drug-repurposing screen (Methods 4.8).

Four modules (DAM-up, IFN-up, chemokine-up, homeostatic-down) queried
against LINCS L1000 chemical-perturbation up/down libraries via the
Enrichr API. Reversal score combines adjusted p-values: disease-up modules
should enrich in a compound's down-regulated genes and vice versa
(weights 1.0/0.8/0.8/1.0).
"""
import json, time
import os
import urllib.request, urllib.parse
import numpy as np
import pandas as pd

OUT = os.environ.get('OUT_DIR', 'out')

MODULES = {
 'DAM_MG': ['APOE','SPP1','LGALS3','CST7','TREM2','TYROBP','GPNMB','ITGAX','CD9','CTSB','CTSD','CSF1'],
 'IFN_module': ['IFIT3','ISG15','USP18','IFIT1','OASL2','IRF7','STAT1'],
 'Chemokine': ['CCL5','CCL3','CCL4','CXCL10','CCL2','CCL12'],
 'Homeostatic_MG': ['P2RY12','TMEM119','CX3CR1','SALL1','SIGLECH','HEXB','FCRLS','OLFML3','GPR34','P2RY13'],
}
# DAM/IFN/Chemokine are disease-UP (want enrichment in compound-DOWN lists);
# Homeostatic is disease-DOWN (want enrichment in compound-UP lists).
DIRECTION = {'DAM_MG': 'down', 'IFN_module': 'down',
             'Chemokine': 'down', 'Homeostatic_MG': 'up'}
WEIGHTS = {'DAM_MG': 1.0, 'IFN_module': 0.8,
           'Chemokine': 0.8, 'Homeostatic_MG': 1.0}

ENRICHR = 'https://maayanlab.cloud/Enrichr'

def enrich(genes, library):
    data = urllib.parse.urlencode(
        {'list': '\n'.join(genes), 'description': 'module'}).encode()
    req = urllib.request.Request(f'{ENRICHR}/addList', data=data)
    uid = json.load(urllib.request.urlopen(req))['userListId']
    url = f'{ENRICHR}/enrich?userListId={uid}&backgroundType={urllib.parse.quote(library)}'
    return pd.DataFrame(json.load(urllib.request.urlopen(url))[library],
                        columns=['rank','term','p','z','combined','genes','q','oldp','oldq'])

LIB_UP = 'LINCS_L1000_Chem_Pert_up'
LIB_DN = 'LINCS_L1000_Chem_Pert_down'

hits = {}
for mod, genes in MODULES.items():
    lib = LIB_DN if DIRECTION[mod] == 'down' else LIB_UP
    df = enrich(genes, lib)
    df['module'] = mod
    hits[mod] = df[['term', 'p', 'q', 'module']]
    time.sleep(1)

allhits = pd.concat(hits.values())
piv = allhits.pivot_table(index='term', columns='module', values='q')
for m in MODULES:
    if m not in piv:
        piv[m] = 1.0
piv = piv.fillna(1.0)
piv['reversal_score'] = -sum(WEIGHTS[m] * np.log10(piv[m].clip(1e-16)) for m in MODULES)
top = piv.sort_values('reversal_score', ascending=False).head(30)
top.to_csv(f'{OUT}/lincs_reversal_top30.csv')
print(top.head(15))
