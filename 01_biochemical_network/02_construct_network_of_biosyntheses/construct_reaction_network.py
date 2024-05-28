import networkx as nx
import pandas as pd
import pickle5 as pickle
import sys
sys.path.append('..')
from utils import *
from datetime import date
today = date.today()
date = today.strftime('%d%b%Y')
print ('Date prefix:', date)

PRECOMPUTED_GRAPH_PATH = ''
REACTION_FILES = ['../parse_reaction_dbs/rhea/08Mar2023_rhea_reaction_smiles_no_cofs.csv', 
                  '../parse_reaction_dbs/bkms/1Sep2023_bkms-mapped.txt',  
                  '../parse_reaction_dbs/metacyc/08Mar2023_metacyc_reaction_smiles_no_cofs.csv']

REACTION_COLUMN = ['reaction_smiles_no_cofs', 'smiles', 'reaction_smiles']

REACTION_SET_NAME = ['_rhea','_bkms','_metacyc']

# Generate network 
try:
    with open(PRECOMPUTED_GRAPH_PATH, 'rb') as f:
        met = pickle.load(f)
    print ('Loading pre-computed graph')
except:    
    print ('RECOMPUTING GRAPH')
    reaction_datasets = []

    for f, c in zip(REACTION_FILES,REACTION_COLUMN):
        df = pd.read_csv(f, sep='\t')
        df['smiles'] = df[c]
        if 'Reaction' in df.columns: 
            reversibles = df[df['Reaction'].map(lambda x: '<=>' in str(x))].copy()
            reversibles.loc[:,'smiles'] = reversibles['smiles'].map(lambda x: flip_reaction(x))
            df = df.append(reversibles).reset_index()
        df = df.dropna(subset='smiles')
        df = df[df['smiles'].map(lambda x: '*' not in x)]
        df['smiles'] = df['smiles'].map(lambda x: standardize_reaction_smiles(x))
        reaction_datasets.append(df)
    
    reaction_df = reaction_datasets[0]
    for i in range(len(reaction_datasets)-1):
        print (REACTION_SET_NAME[i:i+2])
        reaction_df = reaction_df.merge(reaction_datasets[i+1], on='smiles', how='outer',
                                        suffixes=REACTION_SET_NAME[i:i+2])

    for i in ['level_0', 'Unnamed: 0.1', 'Unnamed: 0']:
        if i in reaction_df.columns:
            reaction_df = reaction_df.drop(columns=[i])

    all_smiles = reaction_df['smiles'].values
    all_smiles = list(set([str(r) for r in all_smiles if r]))
    all_smiles = [r for r in all_smiles if len(r.split('>'))==3]
    print ('{} reactions in graph'.format(len(all_smiles)))
    
    print ('Computing graph')
    met = construct_pathway_from_list(all_smiles)

print ('Number of nodes', len(met.nodes))
print ('Number of edges', len(met.edges))
print ('Number of connected components', nx.number_connected_components(nx.Graph(met)))


print ('Number of reaction nodes', len([n for n in met.nodes if '>>' in n]))
print ('Number of chemical nodes', len([n for n in met.nodes if '>>' not in n]))

n_t_s_df = pd.read_csv('../../../spectranalysis/13Jul2023_metabolite_name_to_smiles_df.tsv', sep='\t')

name_to_smiles = dict(zip(n_t_s_df['name'], n_t_s_df['smiles']))
smiles_to_name = dict(zip(n_t_s_df['smiles'], n_t_s_df['name']))



#uncomment to save
reaction_df_path = '{}_all_reaction_from{}.csv'.format(date, ''.join(REACTION_SET_NAME))
reaction_graph_path = '{}_whole_metabolic_network_labeled.pkl'.format(date)

reaction_df.to_csv(reaction_df_path, sep='\t')
nx.write_gpickle(met, reaction_graph_path)

print ('Reactions saved to {}'.format(reaction_df_path))
print ('Network saved to {}'.format(reaction_graph_path))