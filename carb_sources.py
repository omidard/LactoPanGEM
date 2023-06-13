#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 19:17:21 2022

@author: omidard
"""

from dgap import m9,exchange_reactions,essential_substrates, total_dir
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import matplotlib.pyplot as plt
import seaborn as sns







carbs = ['cellb_e',
'fru_e',
'gal_e',
'gam_e',
'glc_D_e',
'malt_e',
'man_e',
'raffin_e',
'rib_D_e',
'sbt_D_e',
'glcn_e',
'lcts_e',
'sucr_e',]



directory = '/home/omidard/gems/allgems/'
models =glob('%s/*.json'%directory)
ex = []
model = load_json_model(models[1])
for reaction in model.reactions:
    if 'EX_' in reaction.id:
        ex.append(reaction.id)
carbs_ex=[]
for i in ex:
    for x in carbs:
        if x in i:
            carbs_ex.append(i)
print(carbs_ex)
growth_rates = pd.DataFrame()
growth_rates['carbon_sources'] = carbs_ex
growth_rates = growth_rates.set_index(growth_rates['carbon_sources'])
growth_rates = growth_rates.drop('carbon_sources', axis=1)
for mod in models:
    model = load_json_model(mod)
    model = m9(model)
    model.reactions.get_by_id('EX_glc_D_e').lower_bound = 0
    gr = []
    for i in carbs_ex:
        try:
            model.reactions.get_by_id(i).lower_bound = -25
            gr.append(model.optimize().fluxes.BIOMASS2)
            model.reactions.get_by_id(i).lower_bound = 0
        except:
            gr.append(0)
    colname=model.id
    name = colname.replace('.json','')
    growth_rates[name]=gr
growth_rates.to_csv('/home/omidard/carbon_sources.csv',index=True)
























