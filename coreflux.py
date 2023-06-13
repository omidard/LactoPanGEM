#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 17:51:48 2022

@author: omidard
"""


from dgap import m9,coreflux
import cobra
import pandas as pd
import multiprocessing as mp
from cobra.io import load_json_model
from glob import glob



models =glob('%s/*.json'%'/home/omidard/gems/allgems/')
pool = mp.Pool(mp.cpu_count())
dfc = pool.map(coreflux, [mod for mod in models])
coreflux = pd.concat(dfc,axis=1)
coreflux.fillna(0, inplace=True)
coreflux.to_csv('/home/omidard/coreflux.csv')


    














