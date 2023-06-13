#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:02:54 2022

@author: omidard
"""


from dgap import m9,exchanges_fluxes,total_dir
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import math




directory = '/home/omidard/gems/allgems/'
models =glob('%s/*.json'%directory)
flux_collector=[]
for mod in models:
    model = load_json_model(mod)
    try:
        m9(model)
        col=model.id
        name = col.replace('.json','')
        df = pd.DataFrame({name:model.optimize().fluxes})
        flux_collector.append(df)
    except:
        pass
species_all_fluxes = pd.concat(flux_collector, axis=1)
species_all_fluxes.to_csv('/home/omidard/species_all_fluxes.csv',index=True)



