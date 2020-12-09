# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 15:32:29 2020

@author: Bátora Dániel
"""
import os 
import numpy as np 
import requests


def get_patent_ids(cid):
    url_start = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    url_end = "/xrefs/PatentID/JSON"
    url = url_start + str(cid) + url_end
    
    r = requests.get(url)
    json_file = r.json()
    
    return json_file["InformationList"]["Information"][0]["PatentID"]





def get_disease(cid, sensitivity = "curated"): 
    
    """
    Comparative Toxicogenomics Database (CTD) based search for 
    disease indications based on PubChem Compound ID (cid)
    
    Parameters: 
        sensitivity: 
            curated: return hits with direct evidence 
            all : return all hits 
    
    Returns: 
        Dictionary of markers and therapeutic indications for the compound
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
                    therapeutic.append(ctd_json[i]["DiseaseName"])
                elif "marker" in ctd_json[i]["DirectEvidence"]: 
                    markers.append(ctd_json[i]["DiseaseName"])
        return {"markers" : markers, "therapeutic": therapeutic}

                
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
    
    
def get_cids_of_disease(disease, db):
    
    """
    PubChem CIDs associated with a given disease 
    
    Parameters: 
        disease(str): name of disease for example "Parkinson"
        db(list): dataset containing compounds with disease indications generated from get_disease
    """
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
    
    """
    Return compounds that have a specific therapeutic
    indication and a specific side effect
    """
    
    thera = get_cids_of_disease(thera, db)
    marker = get_cids_of_marker(marker, db)
    
    return list(set(thera).intersection(marker))