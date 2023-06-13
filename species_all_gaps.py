#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 12:51:12 2022

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



dirs = total_dir('/home/omidard/allgems')
for dec in dirs:
    directory = dec+'/corr'
    genera = dec.replace('/home/omidard/allgems/','')
    models =glob('%s/*.json'%directory)
    gap_collector=[]
    for mod in models:
        gaps=[]
        model = load_json_model(mod)
        for re in model.genes.GAP.reactions:
            gaps.append(re.id)
        df = pd.DataFrame({model.id:gaps})
        gap_collector.append(df)
    species_all_gaps = pd.concat(gap_collector, axis=1)
    species_all_gaps.columns = pd.MultiIndex.from_arrays([[genera for i in range(len(species_all_gaps.columns))], species_all_gaps.columns])
    species_all_gaps.to_csv(dec+'/species_all_gaps.csv')