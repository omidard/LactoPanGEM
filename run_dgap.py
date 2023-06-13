#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 22:22:07 2022

@author: omidard
"""
from dgap import scan,tempfind,gaps,zx,fluxanalyze,addgaps,m9
import multiprocessing as mp
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
import matplotlib.pyplot as plt
from cobra import Model, Reaction, Metabolite
import multiprocessing as mp
import numpy as np


#address to directories
dir1='/home/omidard/allgems/Limosilactobacillus_fermentum/biomassed'
temp_name = 'GCF_009556455.1.json'

output = scan(dir1)
all_gems=output[0]
failed=output[1]
temps=output[2]
temps_inf = tempfind(all_gems,failed,temps)
allmissing = gaps(failed,temp_name)
out = zx(allmissing)
xlist = out[0]
modelsz= out[1]
allgapz=[]
for i in xlist:
    gapz=[]
    pool = mp.Pool(mp.cpu_count())
    target = pool.map(fluxanalyze, [rea for rea in i])
    gapz.append(target)
    allgapz.append(gapz)
    
    
print('this is allgapz[1]',allgapz[1][0])


for i in range(len(modelsz)):
    dir2='/home/omidard/allgems/Limosilactobacillus_fermentum/biomassed/'
    dir3='/home/omidard/allgems/Limosilactobacillus_fermentum/feasible/'
    dir4='/home/omidard/allgems/Limosilactobacillus_fermentum/failed2/'
    model = load_json_model(dir2+modelsz[i])
    template = load_json_model(dir2+temp_name)
    for x in allgapz[i][0]:
        #print('model id  = ',modelsz[i],'\n','gaps = ',allgapz[i][0],'\n','len = ',len(allgapz[i][0]))
        if x != 'no':
            reaction = template.reactions.get_by_id(x)
            reaction2 = Reaction(reaction.id)
            reaction2.name = reaction.name
            reaction2.subsystem = reaction.subsystem
            reaction2.lower_bound = reaction.lower_bound
            reaction2.upper_bound = reaction.upper_bound
            reaction2.add_metabolites(reaction.metabolites)
            reaction2.gene_reaction_rule = '(GAP)'
            model.add_reactions([reaction2])
            model.repair()
            print(reaction2.id)
            print('-----done')
            m9(model)
            print(model.optimize().fluxes.BIOMASS2)
            cobra.io.json.save_json_model(model,dir2+model.id)
            
for i in range(len(modelsz)):
    dir2='/home/omidard/allgems/Limosilactobacillus_fermentum/biomassed/'
    dir3='/home/omidard/allgems/Limosilactobacillus_fermentum/feasible/'
    dir4='/home/omidard/allgems/Limosilactobacillus_fermentum/failed2/'
    model = load_json_model(dir2+modelsz[i])
    m9(model)
    print(model.id,model.optimize().fluxes.BIOMASS2)
    if model.optimize().fluxes.BIOMASS2 >= 0.1:
        cobra.io.json.save_json_model(model,dir3+model.id)
    if model.optimize().fluxes.BIOMASS2 < 0.1:
        cobra.io.json.save_json_model(model,dir4+model.id)
            