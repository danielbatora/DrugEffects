# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:45:42 2020

@author: Bátora Dániel
"""
import pandas as pd
import os 
import pickle

folder = "C:/Users/Bátora Dániel/Downloads/CTD_chemicals.xml"
os.chdir(folder)


import xml.etree.ElementTree as ET

data =  ET.parse("CTD_chemicals.xml")

root = data.getroot()

from tqdm import tqdm 
casr_mesh_dict = {}
for i in tqdm(range(len(root))):
    if root[i][2].tag == "CasRN": 
        casr_mesh_dict.update({root[i][2].text : root[i][1].text.split(":")[-1]})
        
df = pd.DataFrame(casr_mesh_dict.values(), index = casr_mesh_dict.keys())

pickle.dump(casr_mesh_dict, open("casr_mesh_dict.p", "wb"))
df.to_csv("casr_mesh_df.csv")




