# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 14:30:46 2020

@author: Bátora Dániel
"""


import pandas as pd 
import os 
import numpy as np 
import requests
from tqdm import tqdm 
import pickle
folder ="C:/Users/Bátora Dániel/Desktop/Projects/tattoo/Database"

qx_db = pd.read_csv(os.path.join(folder, "quinoxaline.csv"))


syn = qx_db.columns[2]

name = qx_db.columns[1]

syns = [str(i).split("|") for i in qx_db[syn] ] 


names = list(qx_db[name]) 



chloro_qxs = [i for i, n in enumerate(names) if "Chloro" in str(n) or "chloro" in str(n)] 
iodo_qxs = [i for i, n in enumerate(names) if "Iodo" in str(n) or "iodo" in str(n)] 
bromo_qxs = [i for i, n in enumerate(names) if "Bromo" in str(n) or "bromo" in str(n)] 


#search qx_db for disease indications 
hits = []
for i in tqdm(range(len(cids))): 
    if type(cids[i]) == list: 
        for j in cids[i]: 
            try: 
                a = int(ba_db.loc[qx_db["cid"] == j]["cid"])
                hits.append(a)
            except: 
                continue



def casrn_to_cid(casrn):
    url_start = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    url_end = "/cids/JSON"
    
    url = url_start + str(casrn) + url_end
    r = requests.get(url)
    json_file  = r.json()
    try: 
        return json_file["IdentifierList"]["CID"]
    except: 
        return False
        


cids_load = pickle.load(open("cids_with_disease_indications.p", "rb"))

n = len(cids_load)
n=0
cids = cids_load.copy()

for i in tqdm(range(n, len(casr_mesh_dict))): 
    json_file = casrn_to_cid(list(casr_mesh_dict.keys())[i])
    cids.append(json_file)
    pickle.dump(cids, open("cids_with_disease_indications.p", "wb"))

list(casr_mesh_dict.keys())[0]

def get_patent_ids(cid):
    url_start = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    url_end = "/xrefs/PatentID/JSON"
    url = url_start + str(cid) + url_end
    
    r = requests.get(url)
    json_file = r.json()
    
    try: 
        return json_file["InformationList"]["Information"][0]["PatentID"]
    except: 
        return []


def get_disease(cid, sensitivity = "curated"): 
    
    """
    Comparative Toxicogenomics Database (CTD) based search for 
    disease indications based on PubChem Compound ID (cid)
    
    Parameters: 
        sensitivity: 
            curated: return hits with direct evidence 
            all : return all hits 
    
    Returns: 
        List of strings of disease indications 
    """
    
    
    
    assert sensitivity == "curated" or sensitivity == "all"
    url_start = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
    url_end = "/JSON?heading=Associated+Disorders+and+Diseases"
    url = url_start + str(cid)  + url_end
    ctd_url_start = "http://ctdbase.org/tools/batchQuery.go?inputType=chem&inputTerms="

    
    if sensitivity == "curated":
        ctd_url_end = "&report=diseases_curated&format=json"
    else: 
        ctd_url_end = "&report=diseases&format=json"
    
    
    r = requests.get(url)
    json_file = r.json()
    
    try: 
        mesh_ids =  []
        for i in range(len(json_file["Record"]["Reference"])): 
            mesh_id = json_file["Record"]["Reference"][i]["URL"].split("=")[-1]
            mesh_ids.append(mesh_id)
    
    except: 
        return
    
    if sensitivity == "curated":
        markers = []
        therapeutic = []
        for mesh_id in mesh_ids: 
            ctd_url = ctd_url_start + str(mesh_id)  + ctd_url_end
            ctd = requests.get(ctd_url)
            ctd_json = ctd.json()
            for i in range(len(ctd_json)): 
                if ctd_json[i]["DirectEvidence"] == "therapeutic": 
                    therapeutic.append(ctd_json[i]["DiseaseName"].split(":")[-1])
                elif "marker" in ctd_json[i]["DirectEvidence"]: 
                    markers.append(ctd_json[i]["DiseaseName"].split(":")[-1])
        return {"markers" : markers, "therapeutic": therapeutic}, 

                
    elif sensitivity == "all": 
        diseases = []
        
        for mesh_id in mesh_ids: 
           ctd_url = ctd_url_start + str(mesh_id)  + ctd_url_end
           ctd = requests.get(ctd_url)
           ctd_json = ctd.json()
           for i in range(len(ctd_json)): 
               diseases.append(ctd_json[i]["DiseaseName"])
       
        diseases = [str(i) for i in np.unique(np.array(diseases))]
        
        return diseases
    

x = get_disease(2719)
    
    
ba_disease_ids = {}
n = 1
from tqdm import tqdm
for i in tqdm(range(len(hits))):
    try: 
        inds = get_disease(hits[i], sensitivity = "curated")
        if type(inds) == dict: 
            ba_disease_ids.update({hits[i] : inds})
            pickle.dump(ba_disease_ids, open("ba_disease_ids.p", "wb"))
    except: 
        continue




iodo_patents = {}

for i in qx_db["cid"][iodo_qxs]:
    patent_ids = get_patent_ids(i)
    if len(patent_ids) >0:
        print("hit")
        iodo_patents.update({i: patent_ids})




    
    
get_disease(4054, "curated")


def cid_to_casrn(cid): 
    url_start = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
    url_end = "/JSON?heading=Other+Identifiers"
    url = url_start + str(cid)  + url_end
    
    r = requests.get(url)
    try: 
        json = r.json()["Record"]["Section"][0]["Section"][0]["Section"][0]["Information"]
        casr = []
        for i in json: 
            casr.append(i["Value"]["StringWithMarkup"][0]["String"])
        
        
        if all([casr[0] == i for i in casr]): 
            return casr[0]
        else: 
            return casr
    except: 
        return False
i = 0
cas_ids = []
for i in tqdm(range(len(iodo_qxs))): 
    casrn = cid_to_casrn(qx_db["cid"][iodo_qxs[i]])
    if not casrn is False: 
        cas_ids.append(casrn)
    


diseases_iodo_qxs = {}
for i in qx_db["cid"]: 
    hits = get_disease(i)
    if type(hits) == dict: 
        print("HIT")
        diseases_iodo_qxs.update({str(i): hits})
