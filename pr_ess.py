#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 12:27:23 2022
essential genes for certain products
@author: omidard
"""


from dgap import m9,producers,pro_act_ge,pro_ess_ge,pro_ess_eff
import cobra
import pandas as pd
import multiprocessing as mp
from cobra.io import load_json_model
from glob import glob



models =glob('%s/*.json'%'/home/omidard/gems/allgems/')
pool = mp.Pool(mp.cpu_count())
listc = pool.map(producers, [mod for mod in models])
listc2=[]
for l in listc:
    for j in l:
        listc2.append(j)
prod=pd.DataFrame()
prod['id'] = listc2
prod.to_csv('/home/omidard/prod.csv')
        
dfc = pool.map(pro_act_ge, [mod for mod in listc2])
dfc2 = pool.map(pro_ess_ge, [df for df in dfc])
dfc3 = pool.map(pro_ess_eff, [df for df in dfc2])

product_based_single_knock_out = pd.concat(dfc3,axis=1)
pbsko = product_based_single_knock_out.T
wt=[]
for i in pbsko.index:
    model = load_json_model('/home/omidard/gems/allgems/'+i)
    m9(model)
    wt.append(model.optimize().fluxes.EX_mnl_e) #targetmet
pbsko['wild_type'] = wt
pbsko.to_csv('/home/omidard/pbsko.csv')


    
