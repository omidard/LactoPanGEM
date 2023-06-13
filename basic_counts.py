#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 11:08:08 2022

@author: omidard
"""

from dgap import basic_counts
import cobra
import pandas as pd
import multiprocessing as mp
from cobra.io import load_json_model
from glob import glob
import numpy as np


models =glob('%s/*.json'%'/home/omidard/gems/allgems/')
pool = mp.Pool(mp.cpu_count())
dfc = pool.map(basic_counts, [modelid for modelid in models])
counts = pd.concat(dfc,axis=0)
counts.to_csv('/home/omidard/counts.csv')
