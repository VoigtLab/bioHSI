from bs4 import BeautifulSoup
import re
import pandas as pd
import requests
from tqdm import tqdm
import numpy as np
import json

bkms_reactions = pd.read_csv('1Sep2023_bkms-mapped.txt',sep='\t')
kegg_rids = bkms_reactions['Reaction_ID_KEGG'].drop_duplicates().tolist()
outfile = 'kegg_reaction_pages.txt'

def get_reaction_page_for_kegg_rid (krid):
    kegg_gene_url = 'https://www.genome.jp/dbget-bin/www_bget?rn:'
    url = kegg_gene_url+krid
    res = requests.get(url)
    return res.text

unique_kegg_rids =(','.join([str(x) for x in kegg_rids])).split(',')

with open(outfile,'w') as f:
    f.write("")

for krid in tqdm(unique_kegg_rids):
    with open(outfile,'a') as f:
        f.write(get_reaction_page_for_kegg_rid (krid)+"//EOF//")
