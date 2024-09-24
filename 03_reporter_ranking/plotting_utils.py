from datetime import date
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

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
                                           figure_dir=None,
                                           fig_size = (1.5,1.5)):
    
    today = date.today()
    date_str = today.strftime('%d%b%Y')
    subset = uniqueness_scores[uniqueness_scores['{}_steps'.format(organism_abbrev)]<=max_steps].copy()
    subset = subset[subset['lambda_max'] > 200]
    candidates = subset.sort_values(by=metric, ascending=False).head(num_top_cands).loc[:, 
                    ['name', 'smiles','lambda_max','peak_info', metric, '{}_steps'.format(organism_abbrev), 'max_peak']]

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
    to_plot['linewidth'] = np.array([0.4 if to_plot.loc[idx, 'name'] in combined_highlights else 0.1 for idx in to_plot.index]).astype(np.float64)
    to_plot = to_plot.sort_values('size')
    ax = plot_jittery_horizontal_scatterplot(to_plot, x=metric, y='{}_steps'.format(organism_abbrev), size='size', 
                                        sizes={3:primary_size, 2:secondary_size, 1:small_marker_size},
                                        hue='color', palette=palette, 
                                        linewidth=to_plot['linewidth'], edgecolor='black', legend=False)

    ax.set_yticks(list(range(0, max_steps+1, 5))+[max_steps+1])
    ax.set_yticks(range(max_steps+1), minor=True)

    ax.set_title(organism, fontsize = 7)
    ax.set_ylabel('Enzyme steps')
    ax.set_xlabel('Uniqueness')
#     ax.set_yticks(range(12))
    ax.set_yticklabels(list(range(0, max_steps+1, 5))+['>{}'.format(max_steps)])
    ax.tick_params(axis="y", which="minor", right=False, left=True)
    plt.minorticks_on()
    if xlim is not None:
        ax.set_xlim(*xlim)
    
    if figure_dir is not None:
        plt.savefig('{}/{}_steps_vs_{}_uniqueness_{}_v2.pdf'.format(figure_dir, date_str, metric, bb_set_name))