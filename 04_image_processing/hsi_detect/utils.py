from spectral import *
import numpy as np
import spectral.io.envi as envi
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
import cv2 as cv
import scipy
from PIL import Image
from tqdm import tqdm
import PIL
import scipy
from sklearn import mixture
import os 
import re
from collections import OrderedDict
from joblib import Parallel, delayed

from scipy.ndimage import uniform_filter
from scipy.optimize import curve_fit
from matplotlib.patches import ConnectionPatch

from sklearn.cluster import MiniBatchKMeans, AgglomerativeClustering

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def get_closest_wl_ind (centers, wl):
    """
    centers (list of float): 
    wl (float): target wavelength
    
    Returns:
    index of wavelength in centers that is closest to wl
    """
    centers = np.array(centers)
    return np.argmin(np.abs(centers-wl))

def get_ratio (lib, wl_1, wl_2):
    """
    lib: hyperspectral image 
    wl_1 (int): wavelength of numerator
    wl_2 (int): wavelength of denominator
    """
    centers = lib.bands.centers
    ind_1 = get_closest_wl_ind(centers, wl_1)
    ind_2 = get_closest_wl_ind(centers, wl_2)
    
    return lib.read_band(ind_1) / lib.read_band(ind_2)

def get_bands_by_wl (lib, wl_1, wl_2):
    """
    lib: hyperspectral image 
    wl_1 (int): lower bound on wavelength
    wl_2 (int): upper bound on wavelength
    """
    centers = lib.bands.centers
    
    #create HSI with data between designated wavelengths
    wavelength_inds = range (get_closest_wl_ind(centers, wl_1), get_closest_wl_ind(centers, wl_2) + 1)
    
    img = lib.read_bands(wavelength_inds)
    return img

def get_derivative (lib, wl_1, wl_2):
    """
    lib: hyperspectral image 
    wl_1 (int): lower bound on wavelength
    wl_2 (int): upper bound on wavelength
    """
    
    img = get_bands_by_wl (lib, wl_1, wl_2)
    deriv = np.gradient(img, axis=2)
    
    return deriv


def flatten_array (array):
    return np.reshape (array, (array.shape[0]*array.shape[1]))

def cos_sim (x, y):
    """
    x and y are np.arrays of equal length
    """
    
    return np.sum(x * y) / (np.sqrt(np.sum(np.square(x))) * np.sqrt(np.sum(np.square(y))))

# Code for generating reports

def make_mask (pixel_coords, shape):
    mat = np.zeros(shape)
    for pixel in pixel_coords:
        mat[pixel[0], pixel[1]] = 1
    return mat

def select_rectangles (total_shape, list_of_tl_br):
    mask = np.zeros(total_shape) 
    for tl, br in list_of_tl_br:
        mask[tl[0]:br[0], tl[1]:br[1]] = 1
    
    return mask

def hill_eqn(x, ymax, ymin, K, n):
    #internal hill eqn for used for fitting data
    return ymin + (ymax - ymin)*(K**n)/( (x**n) + (K**n))


import inspect


def plot_fit_curve (conc, ab, _hill_eqn, xlims, ylims, figname, markersize=13, fit_linewidth=0.5, xlabel='', 
                    ylabel = 'Fraction absorbed', pts = 5000, ymax_bound = 1, ymin_bound=0,
                    line_color = 'black',
                    marker_color='black', edgecolor='black',
                    width_ratios = [1,10], figsize=(1.5,1.5), ylog=False, error=None,
                    additional_concs_to_plot=None, 
                    additional_data_to_plot = None, 
                    additional_data_ylims = None,
                    colors_for_additional_data=None, x_offset=1e-7, yticks=None, fit_lower_bounds=[-np.inf,-np.inf,0,0], 
                    fit_upper_bounds=[np.inf,np.inf,np.inf, np.inf],
                    ignore_for_fit = None
                   ):
    
    
    f, (ax, ax2) = plt.subplots(1, 2, figsize = figsize, sharey=True, gridspec_kw={'width_ratios': width_ratios}, dpi=300)
    try:
        print ('FITTING')
        
        conc_to_fit = conc
        ab_to_fit = ab
        if ignore_for_fit is not None:
            ignore_idxs=np.ones_like(conc).astype(bool)
            ignore_idxs[ignore_for_fit] = False

            conc_to_fit = conc[ignore_idxs]
            ab_to_fit = ab[ignore_idxs]
        
        params, pcov = curve_fit(_hill_eqn, 
                                             conc_to_fit, ab_to_fit, 
                                             bounds=(fit_lower_bounds, fit_upper_bounds), 
                                             maxfev=10000)
        xs = np.logspace(*np.log10(xlims), pts)#max([int(xlims[1]*10), pts]))
        ydata = _hill_eqn(xs, *params)
        print ('DONE FITTING')
        ax2.plot(xs, ydata, color=line_color, clip_on=False, linewidth=fit_linewidth, zorder=0)
        func_parameters = inspect.signature(_hill_eqn).parameters
        print (['{}: {:.4f}'.format(n.name, x) for x, n in zip(params, list(func_parameters.values())[1:])])
        
        
    except RuntimeError as e:
        print (e)
        pass
    
    conc = np.array(conc)
    
    if error is not None:
        ax.errorbar(np.min(conc), ab[np.argmin(conc)], yerr=[error[np.argmin(conc)]], color='black', linestyle='', linewidth=1, capsize=2)
        ax2.errorbar(conc, ab, yerr=error, capsize=2, color='black', linestyle='', linewidth=1)
    ax.set_xscale('symlog')
    ax2.set_xscale('log')
    
    if ylog:
        ax.set_yscale('symlog')
        ax2.set_yscale('symlog')

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    
    ax2.yaxis.set_visible(False)
    
    linewidth = 0.5
    d = 2.5
    kwargs = dict(marker=[(-1, -d), (1, d)], 
                  markersize=5, 
                  linestyle="none", 
                  color='black', 
                  mec='black', mew=linewidth, clip_on=False)
    
    ax.plot([1.03], [0], transform=ax.transAxes, **kwargs, zorder=100)
    ax2.plot([0], [0], transform=ax2.transAxes, **kwargs, zorder=100)
    

    ax.set_ylabel(ylabel)
    f.text(0.6, 0.01,xlabel, ha='center')
    
    try:
        con = ConnectionPatch(xyA=(0, _hill_eqn(0, *params)), xyB=(xlims[0], np.min(ydata)), coordsA="data", coordsB="data",
                      axesA=ax, axesB=ax2, color=line_color, linestyle='dotted', linewidth=fit_linewidth,
                        connectionstyle='arc3', patchA=None, patchB=None, shrinkA=0.0, 
                        shrinkB=0.0, mutation_scale=1.0, mutation_aspect=None, clip_on=False, zorder=1)
        
        f.add_artist(con)
    except Exception as e:
        print ('Failed', e)
        pass
    
    main_marker_color = marker_color
    if additional_data_to_plot is not None:
        ax3 = ax.twinx()
        ax4 = ax2.twinx()
        for i, (d, ac) in enumerate(zip(additional_data_to_plot, additional_concs_to_plot)):
            if colors_for_additional_data is not None:
                marker_color = colors_for_additional_data[i]
            ax3.scatter(ac[ac==np.min(ac)], d[ac==np.min(ac)], 
                        c = marker_color, s = markersize/2, clip_on=False, 
                        alpha=1, linewidth=markersize/20, zorder=1, marker='s')
            ax3.yaxis.set_visible(False)
            ax3.set_axis_off()
            ax4.sharey(ax3)
            ax4.scatter(ac, d, c = marker_color, s = markersize/2, clip_on=False, 
                        alpha=1, linewidth=markersize/20, zorder=1, marker='s')
            ax4.set_frame_on(False)
            ax4.set_ylim(*additional_data_ylims)
            
    ax.scatter(conc[conc==np.min(conc)], ab[conc==np.min(conc)], c = main_marker_color, 
               s = markersize, clip_on=False, edgecolor=edgecolor, alpha=1, 
               linewidth=markersize/20, zorder=5)
    ax2.scatter(conc, ab, c = main_marker_color, s = markersize, clip_on=False, 
                edgecolor=edgecolor, alpha=1, linewidth=markersize/20, zorder=10)    
    y_offset = 0.0006

    con2 = ConnectionPatch(xyA=(1, 1 - y_offset), xyB=(0, 1 - y_offset), coordsA='axes fraction', coordsB='axes fraction',
                      axesA=ax, axesB=ax2, color="black", linewidth = linewidth,
                        connectionstyle='arc3', patchA=None, patchB=None, shrinkA=0.0, 
                           shrinkB=0.0, mutation_scale=1.0, mutation_aspect=None, clip_on=False)
    f.add_artist(con2)
    
    ax.set_xticks([0])
    
    if yticks is not None:
        ax.set_yticks(yticks)
    ax.set_xticks([], minor=True)
    x_min_exp = int(np.log10(xlims[0]))
    x_max_exp = int(np.log10(xlims[1]))
    step_size = 1
    minor_step_size = 5
    major_ticks  = np.logspace(x_min_exp-1,x_max_exp+1, int((x_max_exp-x_min_exp)/step_size) + 3, base=10)
    minor_ticks = np.concatenate([np.linspace(major_ticks[i], major_ticks[i+1], minor_step_size+1) for i in range(len(major_ticks)-1)])
    ax.minorticks_on()
    ax2.set_xticks(major_ticks)
    ax2.set_xticks(minor_ticks, minor=True,)
    ax2.tick_params(which='minor', length=1.5, width=0.5)
    ax2.tick_params(which='both', width=0.5)
    ax.set_xlim(0-x_offset*markersize, 1e-5)
    ax.set_xticks([], minor=True)
    ax2.set_xlim(*xlims)
    ax2.set_ylim(*ylims)
    
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.08, hspace=None)

    f.savefig(figname, dpi=400)

    plt.show() 


def plot_fit_curves (ls_conc, ls_ab, _hill_eqn, xlims, ylims, figname, markersize=13, fit_linewidth=1, xlabel='Concentration', 
                    ylabel = 'Fraction absorbed', pts = 5000, ymax_bounds = [1], 
                    marker_colors='black', width_ratios = [1,10], figsize=(1.5,1.5), 
                     ylog=False, error=None, x_offset = 3e-6, y_minorticks=True, match_line_colors=True):
    

    
    f, (ax, ax2) = plt.subplots(1, 2, figsize = figsize, sharey=True, gridspec_kw={'width_ratios': width_ratios}, dpi=300)
    linewidth = 0.5
    d = 2.5
    y_offset = 0.0006
    #draw connector at the top of the graph

    kwargs = dict(marker=[(-1, -d), (1, d)], 
                  markersize=5, 
                  linestyle="none", 
                  color= 'black', 
                  mec='black', mew=linewidth, clip_on=False)

    ax.plot([1.03], [0], transform=ax.transAxes, **kwargs, zorder=100)
    ax2.plot([0], [0], transform=ax2.transAxes, **kwargs, zorder=100)
    
    for conc, ab, ymax_bound, marker_color in zip(ls_conc, ls_ab, ymax_bounds, marker_colors):
        try:
            print ('FITTING')
            (ymin, ymax, K, n), pcov = curve_fit(_hill_eqn, conc, ab, bounds=([0,0,0,0], [ymax_bound, ymax_bound, np.inf, np.inf]), maxfev=10000)
            xs = np.logspace(*np.log10(xlims), pts)#max([int(xlims[1]*10), pts]))
            ydata = [_hill_eqn(x, ymin, ymax,  K, n) for x in xs]
            print (f'~~~\nymax:{ymax}\nymin:{ymin}\nK:{K}\nn:{n}\n~~~')
            print ('DONE FITTING')
            line_color = 'black'
            if match_line_colors:
                line_color=marker_color
            ax2.plot(xs, ydata, color = line_color, clip_on=False, linewidth=fit_linewidth, zorder=0)
        except RuntimeError as e:
            print (e)
            pass

        conc = np.array(conc)
        ax.scatter(conc[conc==np.min(conc)], ab[conc==np.min(conc)], c = marker_color, s = markersize, clip_on=False, edgecolor='black', alpha=1, linewidth=markersize/20, zorder=5)
        ax2.scatter(conc, ab, c = marker_color, s = markersize, clip_on=False, edgecolor='black', alpha=1, linewidth=markersize/20, zorder=10)
        ax2.set_xticks([])
        if error is not None:
            ax.errorbar(np.min(conc), ab[np.argmin(conc)], yerr=[error[np.argmin(conc)]], color='black', linestyle='', linewidth=1, capsize=2)
            ax2.errorbar(conc, ab, yerr=error, capsize=2, color='black', linestyle='', linewidth=1)

            
        try:
            #draw connector for line of best fit
            
            ydata_at_min_conc = ydata[np.argmin(np.abs(xs-np.min(conc)))]
            con = ConnectionPatch(xyA=(0, _hill_eqn(0, ymin, ymax,  K, n)), xyB=(xlims[0],ydata_at_min_conc), coordsA="data", coordsB="data",
                          axesA=ax, axesB=ax2, color=line_color, linestyle='dotted', linewidth=fit_linewidth,
                          connectionstyle='arc3', patchA=None, patchB=None, 
                          shrinkA=0.0, shrinkB=0.0, mutation_scale=1.0, 
                          mutation_aspect=None, clip_on=False, zorder=2)
            f.add_artist(con)
        except:
            pass

    ax.set_ylabel(ylabel)
    ax.set_xscale('symlog')
    ax2.set_xscale('log')

    if ylog:
        ax.set_yscale('symlog')
        ax2.set_yscale('symlog')

    
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.yaxis.set_visible(False)

    con2 = ConnectionPatch(xyA=(1, 1 - y_offset), xyB=(0, 1 - y_offset), coordsA='axes fraction', coordsB='axes fraction',
                      axesA=ax, axesB=ax2, color="black", linewidth = linewidth,
                      connectionstyle='arc3', patchA=None, patchB=None, shrinkA=0.0, 
                      shrinkB=0.0, mutation_scale=1.0, mutation_aspect=None, clip_on=False)
    f.add_artist(con2)
    ax.set_xticks([0])
    x_min_exp = int(np.log10(xlims[0]))
    x_max_exp = int(np.log10(xlims[1]))
    step_size = 1
    minor_step_size = 5
    major_ticks  = np.logspace(x_min_exp-1,x_max_exp+1, int((x_max_exp-x_min_exp)/step_size) + 3, base=10)
    minor_ticks = np.concatenate([np.linspace(major_ticks[i], major_ticks[i+1], minor_step_size+1) for i in range(len(major_ticks)-1)])
    print (int((x_max_exp-x_min_exp)/step_size) + 1, major_ticks, minor_ticks)
    if y_minorticks:
        ax.minorticks_on()
    ax2.set_xticks(major_ticks)
    ax2.set_xticks(minor_ticks, minor=True,)
    ax2.tick_params(which='minor', length=1.5, width=0.5)
    ax2.tick_params(axis='both', which='major', width=0.5)
    ax2.tick_params(axis='y', which='both',right=False, labelright=False)
    ax.tick_params(axis='y', which='both',right=False, labelright=False)
    ax.set_xticks([], minor=True)
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.08, hspace=None)    

    ax2.set_xlim(*xlims)
    ax2.set_ylim(*ylims)
    ax.set_ylim(*ylims)
    ax.set_xlim(0-1e-7*markersize-x_offset, 1e-5)    
    
    f.savefig(figname, dpi=400)

    plt.show()
    
def crop_rectangular (im_array, top_left, bottom_right):
    y_max = bottom_right[1]
    y_min = top_left[1]
    x_max = bottom_right[0]
    x_min = top_left[0]
    
    print (y_min, y_max, x_min, x_max)
    return im_array[y_min:y_max, x_min:x_max, :]
    

import colour

def gau(x, mu, sig1, sig2):
    p_1 = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig1, 2.)))
    p_2 = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig2, 2.)))

    return (x<mu)*p_1+(x>=mu)*p_2

    
def reduce_pixel_dimensionality (list_of_pixels, dims=2, random_state=1):
    """
    list_of_pixels (np.array): Nxd array containing N, d-dimensional pixels
    dims (int): number of dimensions to reduce the representation
    
    return:
    Dimensionality-reduced representation of pixels using PCA. 
    NOTE: pixels that contain NaN values will return NaN.
    """
    pixels = []
    pixel_idxs = []
    for i, p in enumerate(list_of_pixels):
            if not any(np.isnan(p)):
                pixels.append(p)
                pixel_idxs.append(i)
    pixel_array = np.stack(pixels)
    pca = PCA(n_components = dims, random_state=random_state, svd_solver='full')
    reduced = pca.fit_transform(pixel_array)
    pca_reduced = reduced
    
    result = np.zeros((len(list_of_pixels), dims)) * np.nan
    result[pixel_idxs] = pca_reduced
    return result


def smooth_spectrum(spectrum, window_size=3):
    """
    Smooth a spectrum using a moving average.

    Parameters:
    - spectrum (numpy array): The input spectrum to be smoothed.
    - window_size (int): The size of the moving average window.

    Returns:
    - smoothed_spectrum (numpy array): The smoothed spectrum.
    """
    if window_size % 2 == 0:
        raise ValueError("Window size must be an odd number.")

    half_window = window_size // 2
    smoothed_spectrum = np.convolve(spectrum, np.ones(window_size)/window_size, mode='valid')
    
    # Extend the smoothed spectrum to match the original length
    pad_before = half_window
    pad_after = len(spectrum) - len(smoothed_spectrum) - pad_before
    smoothed_spectrum = np.pad(smoothed_spectrum, (pad_before, pad_after), mode='edge')
    
    return smoothed_spectrum

def smooth_img_spectrum(im_array, window_size):
    """
    Perform sliding window averaging on axis 2 of a 3D NumPy matrix.

    Parameters:
    - matrix: 3D NumPy array
    - window_size: Size of the sliding window

    Returns:
    - A new matrix with sliding window averages applied along axis 2
    """
    if im_array.ndim != 3:
        raise ValueError("Input matrix must be 3D")

    if window_size <= 0 or window_size > im_array.shape[2]:
        raise ValueError("Invalid window size")
        
    output_matrix = uniform_filter(im_array, size=(1, 1, window_size), mode='nearest')

    return output_matrix

def map_cluster_array(array, mapping_dict):
    # Create the mapping array
    max_key = int(np.max([np.nanmax(array),np.nanmax(list(mapping_dict.keys()))]))
    mapping_array = np.zeros(max_key + 2) * np.nan
    array[np.isnan(array)] = max_key+1
    for key, value in mapping_dict.items():
        k=int(key)
        mapping_array[k] = value
    mapping_array[max_key+1] = np.nan
    
    # Convert the array using the mapping array
    return mapping_array[array.astype(int)]

def filter_out_endmembers_with_reference(endmember_ls, cluster_idx_ls, reference_spec, threshold=1):
    filtered_endmembers_ls = []
    filtered_cluster_idx_ls = []
    for endmembers, clusters in zip(endmember_ls, cluster_idx_ls):
        print (f'Before filtering: {len(endmembers)} endmembers')
        norm_endmembers = endmembers / np.nanmax(endmembers, axis=1, keepdims=True)
        if threshold<1:
            e2e_sim_to_ref = np.zeros((len(norm_endmembers), len(norm_endmembers)))
            for j in range(len(norm_endmembers)):
                for i in range(len(norm_endmembers)):
                    if i != j :
                        e2e_sim_to_ref[j,i] = 1-scipy.spatial.distance.cosine(reference_spec, norm_endmembers[j]-norm_endmembers[i])
                    else:
                        e2e_sim_to_ref[j,i] = 0
            filt_endmembers = np.stack(norm_endmembers)[np.nanmax(e2e_sim_to_ref, axis=0)<=threshold, :]
            #return the indices of the endmembers that pass the filter
            filtered_clusters = clusters.copy()
            for i, passed in enumerate((np.nanmax(e2e_sim_to_ref, axis=0)<=threshold)):
                if not passed:
                    filtered_clusters[filtered_clusters==i] = np.nan
            print (f'After filtering: {len(filt_endmembers)} endmembers')
        else:
            filt_endmembers = endmembers
            filtered_clusters = clusters
            
        filtered_endmembers_ls.append(filt_endmembers)        
        filtered_cluster_idx_ls.append(filtered_clusters)
        
    return filtered_endmembers_ls, filtered_cluster_idx_ls

def cluster_based_extract_endmembers(img, n_clusters, 
                                     output_prefix=None, clustering_method=MiniBatchKMeans, 
                                     reduced_dims=10, return_cluster_idxs=False, 
                                     return_flattened_img=False, norm=True, **kmeans_clustering_kwargs):
    
    img_flattened= np.reshape(img, (img.shape[0]*img.shape[1], img.shape[2]))
    if norm:
        img_flattened = img_flattened / np.nanmax(img_flattened, axis=1, keepdims=True)
    pca_pixels = reduce_pixel_dimensionality(img_flattened, dims=reduced_dims)
    
    pixels = []
    pixel_idxs = []
    for i, p in enumerate(pca_pixels):
            if not any(np.isnan(p)):
                pixels.append(p)
                pixel_idxs.append(i)
    pixel_array = np.stack(pixels)
    
    try:
        iter(n_clusters)
    except:
        n_clusters = [n_clusters]
    
    all_endmembers = []
    all_cluster_idxs = []
    for c in n_clusters:

        kmeans = clustering_method(n_clusters=c, random_state=1, **kmeans_clustering_kwargs)
        cluster_idxs = kmeans.fit_predict(pixel_array)
        
        reconstructed_cluster_idxs = np.zeros(len(img_flattened)) * np.nan
        reconstructed_cluster_idxs[pixel_idxs] = cluster_idxs
        reconstructed_clusters_img = np.reshape(reconstructed_cluster_idxs, img.shape[:2])
        all_cluster_idxs.append(reconstructed_cluster_idxs)
        endmembers = []
        for i in np.unique(cluster_idxs):
            if not np.isnan(i):
                endmembers.append(np.nanmean(img_flattened[reconstructed_cluster_idxs==i,:], axis=0))

        if output_prefix is not None:
            np.save (output_prefix+f'{c}-clusters_LS_cluster_map.npy', reconstructed_clusters_img)
            np.save (output_prefix+f'{c}-clusters_LS_endmembers.npy', endmembers)
        
        all_endmembers.append(endmembers)
    
    if return_flattened_img:
        return all_endmembers, all_cluster_idxs, img_flattened
    if not return_cluster_idxs:
        return all_endmembers
    else:
        return all_endmembers, all_cluster_idxs

        
def kmeans_hierarchical_extract_endmembers(img, output_prefix=None, clustering_method=MiniBatchKMeans,
                                           metric = 'cosine', linkage = 'average', distance_threshold = 0.005,
                                         reduced_dims=3, return_cluster_idxs=False, n_clusters = 1000,
                                           filter_threshold = 0.9, reference_spec=None, norm=True, **kmeans_clustering_kwargs
                                          ):        

    em_ls, clust_ls, img_flattened = cluster_based_extract_endmembers(img, n_clusters, reduced_dims=reduced_dims, return_cluster_idxs=True, 
                                                                      return_flattened_img=True, clustering_method=clustering_method, norm=norm, **kmeans_clustering_kwargs)
    em_ls, clust_ls = filter_out_endmembers_with_reference(em_ls, clust_ls, reference_spec, filter_threshold)
    
    #agglomerative cluster endmembers 
    #set distance cutoff 
    ag = AgglomerativeClustering(n_clusters=None, metric=metric, distance_threshold=distance_threshold, linkage=linkage)

    #map between original (kmeans) clusters to new aggl clusters
    ag_clusters = [ag.fit_predict(em) for em in em_ls]

    mappings = [{x:y for x,y in zip(np.unique(clust_ls), ag_c)} for em, ag_c in zip(em_ls, ag_clusters)]

    relabeled_clusts = [map_cluster_array(clust_ls[i], mappings[i]) for i in range(len(clust_ls))]

    new_endmembers = [[np.mean(img_flattened[ag_clust_labels==n],axis=0) for n in np.unique(ag_clust_labels) if not np.isnan(n)] for ag_clust_labels in relabeled_clusts]        
    
    if return_cluster_idxs:
        return new_endmembers, relabeled_clusts
    else:
        return new_endmembers
    
def UCLS(M, U):
    """
    Copied from https://pysptools.sourceforge.io/_modules/pysptools/abundance_maps/amaps.html#UCLS
    
    Performs unconstrained least squares abundance estimation.

    Parameters:
        M: `numpy array`
            2D data matrix (N x p).

        U: `numpy array`
            2D matrix of endmembers (q x p).

    Returns: `numpy array`
        An abundance maps (N x q).
     """
    Uinv = np.linalg.pinv(U.T)
    return np.dot(Uinv, M[0:,:].T).T


def interpolate_points(point1, point2, total_points):
    # Calculate the step size for each dimension
    if total_points == 1:
        print ('Interpolating across 1 point - returning point 1')
        return [point1]
    else:
        step_x = (point2[0] - point1[0]) / (total_points - 1)
        step_y = (point2[1] - point1[1]) / (total_points - 1)

        # Generate intermediary points
        intermediary_points = []
        for i in range(1, total_points - 1):
            x = point1[0] + i * step_x
            y = point1[1] + i * step_y
            intermediary_points.append((x, y))

        return [point1] + intermediary_points + [point2]
    
def make_rectangle_mask (total_shape, list_of_tl_br):
    mask = np.zeros(total_shape[:2])
    for tl, br in list_of_tl_br:
        mask[tl[0]:br[0]+1, tl[1]:br[1]+1] = 1
    
    return mask

def apply_function_to_rectangles (img, list_of_tl_br, fxn):
    returned = []
    for tl, br in list_of_tl_br:
        selected = img[tl[0]:br[0], tl[1]:br[1]]
        returned.append(fxn(selected))
    return returned

def NNLS(M, U):
    """
    NNLS performs non-negative constrained least squares of each pixel
    in M using the endmember signatures of U.  Non-negative constrained least
    squares with the abundance nonnegative constraint (ANC).
    Utilizes the method of Bro.

    Parameters:
        M: `numpy array`
            2D data matrix (N x p).

        U: `numpy array`
            2D matrix of endmembers (q x p).

    Returns: `numpy array`
        An abundance maps (N x q).

    References:
        Bro R., de Jong S., Journal of Chemometrics, 1997, 11, 393-401.
    """
    import scipy.optimize as opt

    N, p1 = M.shape
    q, p2 = U.shape

    X = np.zeros((N, q), dtype=np.float32)
    MtM = np.dot(U, U.T)
    for n1 in range(N):
        # opt.nnls() return a tuple, the first element is the result
        X[n1] = opt.nnls(MtM, np.dot(U, M[n1]))[0]
    return X

def bandpass_rgb_function (hsi_img, band_centers, coeffs=[1,1,1]):
    #https://photo.stackexchange.com/questions/79534/what-are-the-peak-wavelengths-passing-through-a-bayer-filter
    #loosely based on "Nikon D80 with IR filter removed"
    r = 254*np.mean(hsi_img[:,:, get_closest_wl_ind(600, band_centers):get_closest_wl_ind(690, band_centers)], axis=2)
    g = 254*np.mean(hsi_img[:,:, get_closest_wl_ind(500, band_centers):get_closest_wl_ind(600, band_centers)], axis=2)
    b = 254*np.mean(hsi_img[:,:, get_closest_wl_ind(400, band_centers):get_closest_wl_ind(500, band_centers)], axis=2)
    
    r = r.astype(int)
    g = g.astype(int)
    b = b.astype(int)
    
    rgb = np.stack([r*coeffs[0], g*coeffs[1], b*coeffs[2]], axis=2)
    rgb[rgb>=255] = 254
    return rgb.astype(np.uint8)

def mask_ellipse(array, center, radius_x, radius_y):
    """
    Masks out an ellipse within a NumPy array.
    
    Parameters:
        array (numpy.ndarray): The input array.
        center (tuple): The center coordinates of the ellipse (row, column).
        radius_x (int): The radius of the ellipse along the x-axis.
        radius_y (int): The radius of the ellipse along the y-axis.
    
    Returns:
        numpy.ndarray: The masked array with the ellipse.
    """
    mask = np.zeros_like(array, dtype=bool)
    rows, cols = array.shape[:2]
    cx, cy = center
    
    y, x = np.ogrid[:rows, :cols]
    mask[((x - cx) / radius_x) ** 2 + ((y - cy) / radius_y) ** 2 <= 1] = True
    
    masked_array = array*mask
    return masked_array

def make_circle_mask  (mask_shape, centers, radii):
    mask = np.zeros(mask_shape[:2])
    for c, r in zip(centers, radii):
        mask += mask_ellipse(np.ones_like(mask), c, r, r)
    return mask
    