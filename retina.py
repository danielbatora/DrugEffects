# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:51:51 2020

@author: Bátora Dániel
"""


import pandas as pd
import os 
from tqdm import tqdm
import requests
import numpy as np 
folder = "C:/Users/Bátora Dániel/Desktop/Projects/tattoo/Retina"



df = pd.read_csv(os.path.join(folder, "retina.csv"))



df_ther = df.loc[df["Direct Evidence"] == "therapeutic"]

cids = []
for i in tqdm(range(df_ther.shape[0])):
    if type(df_ther["CAS RN"].iloc[i]) == str:
        cids.append(casrn_to_cid(df_ther["CAS RN"].iloc[i]))
    
    


cids = list(np.array(cids).flatten())


for i in cids: 
    if type(i) == bool:
        cids.remove(i)


cids_flat = [item for sublist in cids for item in sublist]


download_pngs(cids_flat, os.path.join(folder, "images"))

diseases = {}
for i in cids_flat: 
    try: 
        inds = get_disease(i)
        diseases.update({i: inds})
    except:
        continue
    
    


