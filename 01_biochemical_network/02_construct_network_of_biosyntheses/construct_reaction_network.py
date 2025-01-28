import networkx as nx
import pandas as pd
import pickle5 as pickle
from network_utils import *
from datetime import date
import os
today = date.today()
date = today.strftime('%d%b%Y')
print ('Date prefix:', date)

PRECOMPUTED_GRAPH_PATH = ''
REACTION_FILES = ['../00_data/processed/rhea/21Nov2023_rhea_reaction_smiles_no_cofs_with_sequences.csv', 
                  '../00_data/processed/bkms/21Nov2023_bkms-mapped_w_seqs.tsv',  
                  '../00_data/processed/metacyc/21Nov2023_metacyc_reaction_smiles_no_cofs_with_sequences.tsv']

REACTION_COLUMN = ['reaction_smiles_no_cofs', 'smiles', 'reaction_smiles']

REACTION_SET_NAME = ['_rhea','_bkms','_metacyc']
GRAPH_SAVEDIR = '../00_data/processed/graph'

if not os.path.exists(GRAPH_SAVEDIR):
    os.mkdir(GRAPH_SAVEDIR)

# Generate network 
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
    
    #merge sequence columns:
    for idx in reaction_df.index:
        any_seq = [x for x in reaction_df.loc[idx, ['sequence', 'sequence_rhea', 'sequence_bkms']] if not pd.isna(x)] 
        if len(any_seq):
            reaction_df.loc[idx, 'sequence'] = any_seq[0]
    
    # get length of AA sequences    
    seq_lengths = []
    for idx, seq in zip(reaction_df.index, reaction_df['sequence'].values):
        if seq == 'SPONTANEOUS' or reaction_df.loc[idx, 'SPONTANEOUS?']=='T' or \
        reaction_df.loc[idx, 'EC_Number']=='SPONTANEOUS' or 'spon' in str(reaction_df.loc[idx, 'Commentary_MetaCyc']).lower() or \
        'spon' in str(reaction_df.loc[idx, 'Commentary_KEGG']).lower():
            seq_lengths.append(0)
        elif seq is None or pd.isna(seq):
            seq_lengths.append(None)
        else:
            seq_lengths.append(len(seq))
    
    reaction_df['seq_length'] = seq_lengths
    df_subset = reaction_df.loc[:, ['smiles', 'seq_length']]
    
    # remove reactions with no SMILES
    df_subset['smiles'] = [str(r) if r else None for r in df_subset['smiles']]
    df_subset = df_subset.dropna(subset='smiles')
    df_subset = df_subset[df_subset['smiles'].map(lambda x : len(x.split('>'))==3)]
    df_subset = df_subset.groupby('smiles').max()
    
    all_smiles = list([str(r) for r in df_subset.index])
    metadata = [{'aa_seq_len':x} for x in df_subset['seq_length']]
    

    print ('{} reactions in graph'.format(len(all_smiles)))
    
    print ('Computing graph')
    
    met = construct_pathway_from_list(all_smiles, metadata = metadata)

print ('Number of nodes', len(met.nodes))
print ('Number of edges', len(met.edges))
print ('Number of connected components', nx.number_connected_components(nx.Graph(met)))


print ('Number of reaction nodes', len([n for n in met.nodes if '>>' in n]))
print ('Number of chemical nodes', len([n for n in met.nodes if '>>' not in n]))

#uncomment to save
reaction_df_path = f'{GRAPH_SAVEDIR}/{date}_all_reaction_from{"".join(REACTION_SET_NAME)}.csv'
reaction_graph_path = f'{GRAPH_SAVEDIR}/{date}_whole_metabolic_network_labeled.pkl'

reaction_df.to_csv(reaction_df_path, sep='\t')
nx.write_gpickle(met, reaction_graph_path)

print ('Reactions saved to {}'.format(reaction_df_path))
print ('Network saved to {}'.format(reaction_graph_path))