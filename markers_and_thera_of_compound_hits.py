# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:36:00 2020

@author: Bátora Dániel
"""

import os
import pickle
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import shutil
import requests
qx = pickle.load(open("quinone_disease_ids.p", "rb"))

markers = []
thera = []
for i in qx:
    markers.append(qx[i]["markers"])
    thera.append(qx[i]["therapeutic"])
    

markers = [[item for sublist in markers for item in sublist]][0]
thera = [[item for sublist in thera for item in sublist]][0]



marker_unique = np.unique(np.array(markers))
thera_unique = np.unique(np.array(thera))


from collections import Counter



marker_unique = Counter(markers)

df_markers = pd.DataFrame(dict(marker_unique).values(), index = dict(marker_unique).keys())

df_markers_sorted = df_markers.sort_values(["count"], ascending = False)


df_markers.columns = ["count"]

df_markers_sorted.iloc[:15, -1].plot(kind = "bar", color = "k")
plt.ylabel("Count")
plt.title("Quinoline compounds")
plt.tight_layout()
plt.savefig("Markers.png", dpi = 300)

thera_unique = Counter(thera)

df_thera = pd.DataFrame(dict(thera_unique).values(), index = dict(thera_unique).keys())
df_thera.columns = ["count"]

df_thera_sorted = df_thera.sort_values(["count"], ascending = False)



df_thera_sorted.iloc[:15, -1].plot(kind = "bar", color = "k")
plt.ylabel("Count")
plt.title("Quinoline compounds")
plt.savefig("Thera.png", dpi = 300)


def get_cids_of_disease(disease, db):
    ids =[]
    for i in db: 
        if any([disease in j for j in db[i]["therapeutic"]]): 
            ids.append(i)
            
    return ids

def get_cids_of_marker(marker, db):
    ids =[]
    for i in db: 
        if any([marker in j for j in db[i]["markers"]]): 
            ids.append(i)
            
    return ids

seizure = get_cids_of_marker("Seizure", qx)
parkinsons = get_cids_of_disease("Parkinson", qx)
malaria = get_cids_of_disease("Malaria", qx)



def download_pngs(compounds:list, dirname): 
    url_start = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    url_end = "/PNG"
    
    if not os.path.isdir(dirname): 
        os.mkdir(dirname)
    
    for i in compounds: 
        url = url_start + str(i) + url_end
        response = requests.get(url)
        with open(dirname + "/" + str(i) + ".png", 'wb') as f:
            f.write(response.content)
            

def thera_marker_match(db, thera, marker):
    
    thera = get_cids_of_disease(thera, db)
    marker = get_cids_of_marker(marker, db)
    
    return list(set(thera).intersection(marker))


thera_marker_match(qx, "Seizure", "Seizure")

                
            
            
download_pngs(malaria, "Malaria")
download_pngs(seizure, "Seizure")