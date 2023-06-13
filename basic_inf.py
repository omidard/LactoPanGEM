#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 10:36:10 2022

@author: omidard
"""

from dgap import m9,gap_extract,none_ess_gap
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os



names = []
rootdir = '/home/omidard/allgems'
for file in os.listdir(rootdir):
    local = os.path.join(rootdir, file)
    if os.path.isdir(local):
        names.append(local)
        
for name in names:
    directory1 = name+'/corr/'
    models =glob('%s/*.json'%directory1)
    basic_inf = pd.DataFrame()
    modelid=[]
    gr=[]
    rea=[]
    gaps=[]
    genes=[]
    for mod in models:
        model=load_json_model(mod)
        print(directory1,model.id)
        modelid.append(model.id)
        m9(model)
        try:
            gr.append(model.optimize().fluxes.BIOMASS2)
        except AttributeError:
            gr.append('not feasible')
            pass
        rea.append(len(model.reactions))
        gaps.append(len(model.genes.GAP.reactions))
	genes.append(len(model.genes))
    basic_inf['model_id']=modelid
    basic_inf['growth_rate']=gr
    basic_inf['number_of_reactions']=rea
    basic_inf['number_of_gaps']=gaps
    basic_inf['number_of_genes']=genes
    basic_inf.to_csv(name+'/basic_inf.csv')
        
        
        
    
