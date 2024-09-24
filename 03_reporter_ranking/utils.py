import pandas as pd
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, MolStandardize
from tqdm import tqdm
from queue import PriorityQueue
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import itertools
import re
import hashlib
import seaborn as sns
from datetime import date

def neutralize_atoms(smiles):
    #pulled from http://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return Chem.MolToSmiles(mol)


def standardize_reaction_smiles(reaction_smiles):
    try:
        reactants, products = [neutralize_atoms(Chem.MolToSmiles(Chem.MolFromSmiles(s))) for s in reaction_smiles.split('>>')]
        return '>>'.join([reactants, products])
    except:
        return reaction_smiles
    
    
def standardize_smiles(smiles):
    smiles = smiles.replace("'","")
    if ' ' in smiles or 'R' in smiles or '*' in smiles:
        return smiles
    try:
        smiles = neutralize_atoms(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
        return smiles
    except:
        if Chem.MolFromSmiles(smiles)==None:
            m = Chem.MolFromSmiles(smiles, sanitize=False)
            fix_c = AllChem.ReactionFromSmarts('[#6-:1]>>[C;+0:1]')
            if m:
                f = fix_c.RunReactants([m])
                if len (f)>0:
                    f = f[0][0]
                    print ('fixed smiles :', smiles,Chem.MolToSmiles(f))
                    return neutralize_atoms(Chem.MolToSmiles(f))
                else:
                    return smiles
            else:
                return smiles
        else:
            return smiles
        
def get_peaks(spectrum, wl_index, **kwargs):
    peak_idxs = signal.find_peaks(spectrum, **kwargs)
    return wl_index[peak_idxs[0]]

def metabolic_dijkstra (g, starting_nodes, path_length_field="path_length", 
                        shortest_path_field="shortest_pathway", 
                        ignore_spontaneous=True, aa_seq_len_field='aa_seq_len'):
    """
    Implementation based on Gallo et al, 1993
    """
    g = g.copy()
    pq = PriorityQueue()  
    for node in g.nodes:
        g.nodes[node][path_length_field] = np.inf
        g.nodes[node][shortest_path_field] = []
        g.nodes[node]['in_queue'] = False
        if '>>' in node:
            g.nodes[node]['visit_counter'] = 0
            g.nodes[node]['pred_count'] = len(list(g.predecessors(node)))

    for node in starting_nodes:
        g.nodes[node][path_length_field] = 0
        if '>>' not in node:
            pq.put((g.nodes[node][path_length_field], node))
            g.nodes[node]['in_queue'] = True

    
    while not pq.empty():
        _, curr_node = pq.get() # select nearest node
        g.nodes[curr_node]['in_queue'] = False # remove from queue
        for reaction_node in list(g.successors(curr_node)): # iterate through reactions that use curr node as reactant
            g.nodes[reaction_node]['visit_counter'] += 1 # increment visit counter for reaction node  
            if g.nodes[reaction_node]['visit_counter'] == g.nodes[reaction_node]['pred_count']: # if reaction node has been visited as many times as it has reactants
                f = sum([g.nodes[r][path_length_field] for r in g.predecessors(reaction_node)])
                f_path = list(itertools.chain(*[g.nodes[r][shortest_path_field] for r in g.predecessors(reaction_node)]))
                g.nodes[reaction_node][path_length_field] = f
                g.nodes[reaction_node][shortest_path_field] = f_path
                for p in g.successors(reaction_node): # add all product nodes to queue
                    orig_path_length = g.nodes[p][path_length_field]
                    if  orig_path_length > f + 1: # currently all edge weights = 1
                        if not ignore_spontaneous or (ignore_spontaneous and  g.nodes[reaction_node][aa_seq_len_field]>0):
                            g.nodes[p][path_length_field] = f + 1
                        else: 
                            g.nodes[p][path_length_field] = f 
                        g.nodes[p][shortest_path_field] = f_path + [reaction_node]
                        if not g.nodes[p]['in_queue']:
                            g.nodes[p]['in_queue'] = True
                            pq.put((g.nodes[p][path_length_field], p))
                            if orig_path_length < np.inf:
                                for e in list(g.successors(p)):
                                    g.nodes[e]['visit_counter'] -= 1

    return g


def flip_reaction (reaction_smiles):
    reactants, products = reaction_smiles.split('>>')
    return '>>'.join([products, reactants])

def show_rxn_list (list_of_reaction_smiles, reaction_info_df=[]):
    rd_objs = [AllChem.ReactionFromSmarts(r, useSmiles=True) for r in list_of_reaction_smiles]
    
    
    
    for string, obj in zip(list_of_reaction_smiles, rd_objs):
        plt.figure(figsize=(8,10))
        print(string)
        if len(reaction_info_df) > 0:
            reaction_info = reaction_info_df[reaction_info_df['smiles'] == string].astype(str)
            #stringify
            print (dict(reaction_info.apply(lambda x : ';'.join(x), axis=0)))
        plt.imshow(Chem.Draw.ReactionToImage(obj, subImgSize=(300,500)))
        plt.xticks([])
        plt.yticks([])
        plt.show()


def plot_vs_steps(df, x_col, y_col, figname, xlabel=None, ylabel=None, size_col=None, show=True, color='gray', figsize=(5,5), **kwargs):
    markersize=13
    fit_linewidth=1

    width_ratios = [1,10]

    
    f, (ax, ax2) = plt.subplots(2, 1, figsize = figsize, gridspec_kw={'height_ratios': width_ratios}, sharex=True)

    if size_col:
        ax2.scatter(df[x_col], df[y_col], color = color, edgecolor='black',  zorder=0, s=df[size_col], **kwargs)
    else:
        ax2.scatter(df[x_col], df[y_col], color = color, edgecolor='black',  zorder=0, **kwargs)

    no_known_path = df[df[y_col]==np.inf]
    
    if size_col:
        ax.scatter(no_known_path[x_col], [15]*len(no_known_path), edgecolor='black', zorder=5, color = color, s = no_known_path[size_col], **kwargs)
    else:
        ax.scatter(no_known_path[x_col], [15]*len(no_known_path), edgecolor='black', zorder=5, color = color, **kwargs)

    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    # ax2.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)

    # ax.set_xlim(0-1e-7*markersize, 1e-5)

    linewidth = 0.8
    y_offset = 0#0.0006
    d = -1
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10+linewidth, linestyle="none", color='black', mec='black', mew=linewidth, clip_on=False)

    ax.plot([0], [0], transform=ax.transAxes, **kwargs)
    ax2.plot([0], [1], transform=ax2.transAxes, **kwargs)

    # ax.set_ylabel(ylabel)
    # f.text(0.6, 0.01, 'Concentration ({})'.format(unit), ha='center')
    con = ConnectionPatch(xyA=(1, 0 - y_offset), xyB=(1, 1 - y_offset), coordsA='axes fraction', coordsB='axes fraction',
                      axesA=ax, axesB=ax2, color="black", linewidth = linewidth, 
                        connectionstyle='arc3', patchA=None, patchB=None, shrinkA=0.0, shrinkB=0.0, mutation_scale=1.0, mutation_aspect=None, clip_on=False)
    f.add_artist(con)



    con2 = ConnectionPatch(xyA=(0, 0 - y_offset), xyB=(0, 1 - y_offset), coordsA='axes fraction', coordsB='axes fraction',
                      axesA=ax, axesB=ax2, color="black", linewidth = linewidth,linestyle='dotted',
                        connectionstyle='arc3', patchA=None, patchB=None, shrinkA=0.0, shrinkB=0.0, mutation_scale=1.0, mutation_aspect=None, clip_on=False)
    f.add_artist(con2)


    #     plt.tight_layout(pad=2)

    yticklabels = range(0, int(np.max([x for x in df[y_col].unique() if x < np.inf]))+1, 2)

    ax.set_yticks([15])
    ax.set_yticklabels(['Unknown'])

    ax2.set_yticks(yticklabels)
    ax2.set_yticklabels(yticklabels)
    
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    
    
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.15)

    f.savefig(figname, dpi=400)
    if show:
        plt.show()
    
    return f, ax, ax2

def clean_name(string: str) -> str:
    new = re.sub("[\s\.,]", "_", string)
    new = re.sub("[\[\]\(\)']", "", new)
    new = re.sub('&rarr;', '->', new)
    new = re.sub('<[a-zA-z]*>|<\/[a-zA-z]*>|;|&|^[Aa]n |^[0-9]* ','', new)
    new = re.sub('\+', 'plus', new)
    new = re.sub('^-', 'minus', new)
    new = re.sub(',|\/|\\\\', '-', new)
    return new

def get_smiles_from_name (n, n_t_s_df, g):
    matched = n_t_s_df[n_t_s_df['name']==n].loc[:,'smiles'].drop_duplicates()
    if len(matched)==0:
        matched = n_t_s_df[n_t_s_df['cleaned_name']==clean_name(n)].loc[:,'smiles'].drop_duplicates()    
    return matched.tolist()

def get_name_from_smiles (s, n_t_s_df):
    matched = n_t_s_df[n_t_s_df['smiles']==s].loc[:,'cleaned_name'].drop_duplicates()
    return matched.tolist()

def hash_smiles(smiles):
    """
    Generate modified md5 hash of input SMILES string.
    """
    return str(int(hashlib.md5(smiles.encode("utf-8")).hexdigest(), 16))

def get_canon_tautomer(smiles):
    try:
        std_mol = MolStandardize.rdMolStandardize.CanonicalTautomer(Chem.MolFromSmiles(smiles))
    except:
        std_mol = None
    return std_mol

def plot_jittery_horizontal_scatterplot(data, x, y, jitter=0.25, **kwargs):
    #apply jitter along y axis:
    data = data.copy()
    if 'jitter_multiplier' in data.columns:
        data[y] = data[y] + np.random.uniform(-jitter,jitter,len(data))*data['jitter_multiplier']
    else:
        data[y] = data[y] + np.random.uniform(-jitter,jitter,len(data))
    ax = sns.scatterplot(data, x=x, y=y, **kwargs)
    ax.minorticks_on()
    ax.tick_params(axis='y', which='minor', left=False)
    return ax

def plot_steps_vs_uniqueness_for_organism (uniqueness_scores, bb_set_name, 
                                           organism, organism_abbrev, 
                                           primary_names_to_highlight,
                                           max_steps = 3,
                                           metric = 'global_50',
                                           secondary_names_to_highlight=[],
                                           primary_color='cornflowerblue',
                                           secondary_color='salmon', secondary_size=12,
                                           num_top_cands=10, xlim=None,
                                           primary_size=12, small_marker_size=0.6, 
                                           figure_dir="../../Manuscripts/HSI_manuscript/figure_imgs",
                                           fig_size = (1.5,1.5)):
    
    today = date.today()
    date_str = today.strftime('%d%b%Y')
    subset = uniqueness_scores[uniqueness_scores['{}_steps'.format(organism_abbrev)]<=max_steps].copy()
    subset = subset[subset['lambda_max'] > 200]
    candidates = subset.sort_values(by=metric, ascending=False).head(num_top_cands).loc[:, 
                    ['name', 'smiles','lambda_max','peak_info', metric, '{}_steps'.format(organism_abbrev)]]

    display(candidates)
    
    max_steps = 10
    to_plot = uniqueness_scores.copy()
    to_plot = to_plot[to_plot['lambda_max']>200]
    to_plot.loc[to_plot['{}_steps'.format(organism_abbrev)]>max_steps, '{}_steps'.format(organism_abbrev)] = max_steps+1

    plt.figure(figsize=fig_size, dpi=250)
    to_plot['color'] = [to_plot.loc[idx, 'name'] if to_plot.loc[idx, 'name'] in primary_names_to_highlight else 0 for idx in to_plot.index]
    palette={0: 'lightgray'}
    palette.update({n:primary_color for n in primary_names_to_highlight})
    
    if isinstance(secondary_color,str):
        to_plot['color'] = [to_plot.loc[idx, 'color'] if to_plot.loc[idx, 'name'] not in secondary_names_to_highlight else 2 for idx in to_plot.index]
        palette.update({n:secondary_color for n in secondary_names_to_highlight})
    
    elif isinstance(secondary_color,dict):
        to_plot['color'] = [to_plot.loc[idx, 'color'] if to_plot.loc[idx, 'name'] not in secondary_names_to_highlight else to_plot.loc[idx,'name'] for idx in to_plot.index]
        palette.update(secondary_color)
    
    else:
        raise ValueError('secondary color must be a string or dictionary')
    
    combined_highlights = primary_names_to_highlight + secondary_names_to_highlight
    to_plot['size'] = [3 if to_plot.loc[idx, 'name'] in primary_names_to_highlight else 1 for idx in to_plot.index]
    to_plot['size'] = [2 if to_plot.loc[idx, 'name'] in secondary_names_to_highlight else to_plot.loc[idx,'size'] for idx in to_plot.index]
    to_plot['jitter_multiplier'] = [0.8 if to_plot.loc[idx, 'name'] in combined_highlights else 1 for idx in to_plot.index]
    to_plot['linewidth'] = np.array([0.2 if to_plot.loc[idx, 'name'] in combined_highlights else 0.1 for idx in to_plot.index]).astype(np.float64)
    to_plot = to_plot.sort_values('size')
    ax = plot_jittery_horizontal_scatterplot(to_plot, x=metric, y='{}_steps'.format(organism_abbrev), size='size', 
                                        sizes={3:primary_size, 2:secondary_size, 1:small_marker_size},
                                        hue='color', palette=palette, 
                                        linewidth=to_plot['linewidth'], edgecolor='black', legend=False)
    ax.set_title(organism, fontsize = 7)
    ax.set_ylabel('Enzyme steps')
    ax.set_xlabel('Uniqueness')
    ax.set_yticks(range(12))
    ax.set_yticklabels(list(range(max_steps+1))+['>{}'.format(max_steps)])
    if xlim is not None:
        ax.set_xlim(*xlim)
    plt.savefig('{}/{}_steps_vs_{}_uniqueness_{}_v2.pdf'.format(figure_dir, date_str, metric, bb_set_name))


def min_max_norm(a):
    a = np.array(a).astype(np.float64)
    return (a - np.nanmin(a))/(np.nanmax(a)-np.nanmin(a))