# NOMBRE (Cerda ,  Jean-Pierre): 

from Bio import Entrez
from Bio import SeqIO
from Bio import GenBank
import re
import pandas as pd
import matplotlib.pyplot as plt
import csv as csv
    
def download_pubmed (keyword):
    """
    El siguiente comando permite buscar articulos en pubmed mediante palabras claves
    """
     
    Entrez.email = 'jean.cerda@est.ikiam.edu.ec'
    busq = Entrez.read(Entrez.esearch(db="pubmed", 
                            term=keyword,
                            usehistory="y"))
    webenv = busq["WebEnv"]
    query_key = busq["QueryKey"]
    handle = Entrez.efetch(db="pubmed",
                           rettype="medline", 
                           retmode="text", 
                           retstart=0,
                           retmax=543, webenv=webenv, query_key=query_key)
    data = handle.read()
    dataexp = re.sub(r'\n\s{6}','', data)
    return dataexp
    

