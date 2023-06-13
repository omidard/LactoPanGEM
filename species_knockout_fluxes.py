#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 13:18:49 2022
knock_outs effects
@author: omidard
"""

from dgap import m9,exchanges_fluxes,total_dir, knockout_fluxes
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



dirs = total_dir('/home/omidard/allgems')
for dec in dirs:
    directory = dec+'/corr'
    models =glob('%s/*.json'%directory)
    flux_collector=[]
    for mod in models:
        model = load_json_model(mod)
        for re in model.reactions:
            model = load_json_model(mod)
            model.reactions.get_by_id(re.id).lower_bound = 0
            model.reactions.get_by_id(re.id).upper_bound = 0
        try:
            flux_frame = knockout_fluxes(model)
            flux_frame.columns = pd.MultiIndex.from_arrays([[re.id for i in range(len(flux_frame.columns))], flux_frame.columns])
            flux_collector.append(flux_frame)
        except:
            pass
    species_knockout_exchange_fluxes = pd.concat(flux_collector, axis=1)
    for i in species_knockout_exchange_fluxes.index:
        rows = species_knockout_exchange_fluxes.loc[i]
        if all([ v == 0 for v in rows ]):
            species_knockout_exchange_fluxes = species_knockout_exchange_fluxes.drop(index=(i))
            species_knockout_exchange_fluxes.to_csv(dec+'/species_knockout_exchange_fluxes.csv')
    
    #plot
    big=list(range(0,int(round(math.sqrt(len(species_knockout_exchange_fluxes.index))))+1))
    lst = [list(z) for z in itertools.product(big, repeat=2)]
    fig, axes = plt.subplots(round(math.sqrt(len(species_knockout_exchange_fluxes.index)))+1, round(math.sqrt(len(species_knockout_exchange_fluxes.index)))+1, figsize=(20, 30))
    for i in range(len(species_knockout_exchange_fluxes.index)):
        rows = species_knockout_exchange_fluxes.loc[species_knockout_exchange_fluxes.index[i]].T
        data=pd.DataFrame()
        data[rows.name]=[v for v in rows]
        sns.swarmplot(ax=axes[lst[i][0],lst[i][1]],data=data,palette = "Wistia_r",size=2.7)
        for v in rows:
            if v>0:
                palette = sns.color_palette("magma")
            else:
                palette = "Set1"
        sns.boxplot(ax=axes[lst[i][0],lst[i][1]], data=data, palette=palette)
        fig.savefig(dec+'/knock_out.png',bbox_inches='tight')     
        
        
        
        
        
        
        
        
        
        
        