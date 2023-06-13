#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 11:17:04 2022

@author: omidard
"""

import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
import matplotlib as plt
import seaborn as sns


from dgap import m9,exchange_reactions,essential_substrates, total_dir
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import multiprocessing as mp


cs=pd.read_csv('/home/omidard/minimal_media.csv')
"""
for i in range(len(cs['carbon_sources'])):
    cs['carbon_sources'][i]='growth_on_'+ cs['carbon_sources'][i]
"""
for i in range(len(cs['metabolites'])):
    cs['metabolites'][i]='eliminated_'+ cs['metabolites'][i]
#cs.set_index('carbon_sources',inplace=True)
cs.set_index('metabolites',inplace=True)
cso=cs.T
all_inf=pd.read_csv('/home/omidard/carbon_sources.csv')
#all_inf.drop('Unnamed: 0',axis=1,inplace=True)
#all_inf.drop('Unnamed: 0.1',axis=1,inplace=True)
all_inf.set_index('Unnamed: 0',inplace=True)
carbon_sources = pd.concat([all_inf,cso],axis=1)
carbon_sources.to_csv('/home/omidard/inf.csv',index=True)





