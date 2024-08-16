from bs4 import BeautifulSoup
import re
import pandas as pd
import requests
from tqdm import tqdm
import numpy as np
import json
import pdb
import os

# with open('kegg_reaction_pages.txt','r') as f:
#   kegg_reaction_pages = f.read().split("//EOF//")
outfile = 'kegg_compound_pages.txt'
all_cpds_page = requests.get('https://www.genome.jp/dbget-bin/www_bfind_sub?dbkey=compound&keywords=C&mode=bfind&max_hit=nolimit')
soup = BeautifulSoup(all_cpds_page.text, 'html.parser')
a = soup.find_all('a')
cpd_links = [x.get('href') for x in a if x.get('href')]
if not os.path.isfile(outfile):
  with open(outfile,'w') as f:
    f.write('')

else:
  with open(outfile,'r') as f:
    results = f.read().split('//EOF//')
  
  completed = [x.split('\n')[0] for x in results]


for c in tqdm(cpd_links):
  if c not in completed:
    mol_texts = requests.get('https://www.genome.jp/entry/-f+m+'+c.split('/')[-1]).text
    with open(outfile,'a') as f:
      f.write(c+'\n'+mol_texts+'//EOF//')


# class KeggReactionPage():
   
#     def __init__(self, reaction_page):
#         self.reaction_page = reaction_page
#         self.soup = BeautifulSoup(self.reaction_page, 'html.parser')
#         self.reaction_info_dict = {}
#         self.get_reaction_info()
#         # self.get_substrate_table_info()

#     def get_reaction_info(self):
#         # return a dictionary of info about the reaction
#         header_divs = self.soup.find_all('a')
#         pdb.set_trace()
#         for tab_div in header_divs:
#            if tab_div:
#             for r_div in tab_div.find_all('tr'):
#               if r_div:
#                 if 'Equation' in r_div.text:
#                    self.reaction_info_dict['Equation']  = r_div.text.split('\n')[1]
#                 if 'Other DBs' in r_div.text:
#                    print (r_div, r_div.text)
#                    self.reaction_info_dict['Other DBs']  = r_div.text.split('\n')
#         pdb.set_trace()
                  
            # tab_num = int(re.findall('(?<=tab)[0-9]+(?=_head)', tab_div.get('id'))[0])
            # self.table_headers_dict[tab_num] = [x.strip() for x in tab_div.text.split('\n') if len(x.strip())]

    # def _get_substrates_table_id(self):
    #     #return id for SUBSTRATES table
    #     titles = ['SUBSTRATE','PRODUCT','REACTION DIAGRAM','ORGANISM','UNIPROT']
    #     substrate_table_id = [x for x in self.table_headers_dict if self.table_headers_dict[x][:5]==titles]
    #     return substrate_table_id
    
    # def get_substrate_table_info(self):
    #     substrate_table_id = self._get_substrates_table_id()[0]
    #     hidden_divs = self.soup.find_all('div', id=re.compile('^tab{}'.format(substrate_table_id)))
    #     res_list = []
    #     for div in hidden_divs:
    #         # Extract values from the first, second, and fifth columns
    #         res_dict = {}
    #         for i,c in enumerate(self.table_headers_dict[substrate_table_id]):
    #             cell = div.find('div', class_='cell', id=lambda x: x and x.endswith('c{}'.format(i)))
    #             if cell:
    #                 value = cell.text.strip()
    #                 res_dict[c] = value
    #         if res_dict.keys():
    #             res_list.append(res_dict)
    #     self.substrate_table_info = res_list




# def extract_compound_ids(kegg_reaction_page):
#   soup = BeautifulSoup(kegg_reaction_page, 'html.parser')
#   rid = re.findall('R[0-9]+', soup.find('title').text)[0]
  
#     # Find the tr element with the specific class
#   td = soup.find_all('td')
#   for t in td:
#     if t.get("class"):
#       print (t)
#   # Find all a elements within the tr element
#   a_elements = tr.find_all('a')
  
#   # Extract the identifiers that start with 'C'
#   c_ids = [re.findall('C[0-9]+', a.text)[0] for a in a_elements if re.match('C[0-9]+', a.text)]
  
#   print(rid)
#   print(c_ids)
  

# # extract_compound_ids(kegg_reaction_pages[0])
  
# kr = KeggReactionPage(kegg_reaction_pages[0])