import pandas as pd

from scipy import special, stats
from scipy.spatial.distance import cdist, cosine
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import itertools
from tqdm import tqdm
matplotlib.rcParams['pdf.fonttype'] = 42

# from spectranalysis.utils import *
from joblib import Parallel, delayed

def remove_consecutive_ints (int_ls, method='keep_first'):
    """
    Args:
        int_ls (list[int]): list of ascending integer values
        method (str): `keep_first` or `midpoint`. If `keep_firt`
            keep the lowest value integer in the list of consecutive
            values. If `midpoint`, keep the midpoint value in the list 
            of consecutive integers.
    Returns:
        a list of integers without consecutive values
    """
    res_ls = []
    if method=='keep_first':
        last_number = np.inf
        for i in int_ls:
            if (i - last_number) == 1:
                pass
            else:
                res_ls.append(i)

            last_number = i
    
    elif method=='midpoint':
        last_number = np.inf
        consecutive_lists = []
        for i in int_ls:
            if (i - last_number) == 1:
                consecutive_lists[-1].append(i)
            else:
                consecutive_lists.append([i])
            last_number = i
        res_ls = [ls[int(len(ls)/2)] for ls in consecutive_lists] # find midpoint
    else:
        res_ls = None
        print ('No such method: {}!'.format(method))
    return res_ls

def parse_spectrum (spectrum, intensity_threshold=0.05, derivative_threshold=0.0001, method='midpoint',
                   threshold_subspectra = True, validate=False):
    """
    Split a spectrum at valley points into a list of spectra 
    Args:
        spectrum (array-like): spectrum to parse
        intensity_threshold (float): relative intensity value under which a value in the spectrum
            is considered a valley (separating peaks)
        derivative_threshold (float): value for the derivative of a point with a positive second 
            derivative under which the point is considered to be a valley
        method (string): method to use for picking which point in a valley to split the spectra on
        threshold_subspectra (bool): if true, remove subspectra that do not have any points above the 
            `intensity_threshold`
        validate(bool): if true, raise an error if there is a discrepancy greater than 5% between
            the sum of the area under the subspectra and the area under the original spectrum
    Returns:
        A list of subspectra that sum to `spectrum`. Each subspectrum contains a single 
        identified peak
    """
    
    spectrum = np.array(spectrum)
    norm_spectrum = spectrum / np.max(spectrum)
    ERROR_MARGIN = 0.05*np.sum(norm_spectrum)

    dx = np.gradient(norm_spectrum)
    ddx = np.gradient(dx)
    
    min_idxs = np.where((norm_spectrum<intensity_threshold) | (np.abs(dx) <= derivative_threshold) * (ddx > 0))[0]
    split_points = remove_consecutive_ints(min_idxs, method)
    if len(split_points):
        split_points = [0] + split_points + [-1]
        basis_spectra = []
        
        for i in range(len(split_points)-1):
            partial_spec = spectrum.copy()
            partial_spec[:split_points[i]] = 0
            partial_spec[split_points[i+1]:] = 0
            
            # don't include fake spectra 
            if threshold_subspectra and any(norm_spectrum[split_points[i]:split_points[i+1]]>intensity_threshold):
                basis_spectra.append(partial_spec)
            elif not threshold_subspectra:
                basis_spectra.append(partial_spec)
                
    else:
        return [spectrum]
    if validate:
        assert 1-np.sum(np.abs(np.sum(basis_spectra, axis=0)/np.sum(spectrum))) < ERROR_MARGIN
    return basis_spectra
    
    
def compute_scores (distances, original_spectra, original_to_parsed_idxs_dict):
    """
    returns a uniqueness score for a molecule using a mapping from molecules to peak indexes and a 
    matrix of distances between peaks and background spectra. The score is the sum of the mean of the distances
    between each of a moelcule's peaks and all background peaks.
    """
    scores = []
    for i in range(len(original_spectra)):
        idxs = original_to_parsed_idxs_dict[i]
        score = np.sum(np.nanmean(distances[:,idxs], axis=0))
        scores.append(score)
    return scores

def get_mean_kl_div (spectra, backgrounds, normalize=True):
    if normalize:
        backgrounds = [np.array(b) / np.max(b) for b in backgrounds]
        spectra = [np.array(s) / np.max(s) for s in spectra]
    mean_kl_divs = []
    for s in tqdm(spectra):
        s_kl_divs = []
        for b in backgrounds:
            kl_div = np.sum(special.kl_div(s, b))
            s_kl_divs.append(kl_div)
        mean_kl_divs.append(np.mean(s_kl_divs))
    return mean_kl_divs

def compute_scores_from_subset (distances, original_spectra, original_to_parsed_idxs_dict, fraction_size=0.1):
    """
    returns a uniqueness score for a molecule using a mapping from molecules to peak indexes and a 
    matrix of distances between peaks and background spectra. Returns the sum of the mean of the distances of 
    each peak in a spectrum and all ofn the background spectra. 
    """
    scores = []
    k = int(np.ceil(fraction_size*distances.shape[0]))
    print("Comparing to {} molecules".format(k))
    for i in range(len(original_spectra)):
        idxs = original_to_parsed_idxs_dict[i]
        most_similar_idxs = np.argsort(np.sum(distances[:, idxs], axis=1))[1:1+k]
        if len(idxs) > 0:
            most_similar_idxs = np.repeat(most_similar_idxs[:, np.newaxis], len(idxs), axis=1)
        score = np.sum(np.nanmean(distances[most_similar_idxs, idxs], axis=0))
        scores.append(score)
    return scores

def compute_scores_from_subset_no_parsing (distances, fraction_size=0.1):
    scores = []
    k = int(np.ceil(fraction_size*distances.shape[0]))
    print("Comparing to {} molecules".format(k))
    most_similar_idxs_mat = []
    for i in range(len(distances)):
        # start at index 1 to skip comparison to self
        most_similar_idxs = np.argsort(distances[i,:])[1:1+k] 
        score = np.nanmean(distances[i, most_similar_idxs], axis=0)
        scores.append(score)
        most_similar_idxs_mat.append(most_similar_idxs)
    return np.array(scores), np.stack(most_similar_idxs_mat)

def get_impact_score_ (spectrum, bg_spectra, metric=cosine, intensity_factor=0.1):
    if np.sum(spectrum)==0:
        return 0
    else:
        dists = [metric(bg, bg-spectrum*intensity_factor*bg) for bg in bg_spectra]
        return np.nanmean(dists)

def get_impact_scores(spectra_to_score, bg_spectra, metric=cosine, intensity_factor=0.1, n_jobs=1):
    """
    Score molecule spectral impact on a set of background spectra from a hyperspectral image
    """
    impact_scores = Parallel(n_jobs=n_jobs)(delayed(get_impact_score_)(s, bg_spectra, metric=metric, intensity_factor=intensity_factor) for s in tqdm(spectra_to_score))
    return impact_scores

def get_background_specs(spec_df, metabolites_file, name_to_smiles_dict):
    """
    Get the spectra for the genomic background metabolites 
    """
    background_molecules_df = pd.read_csv(metabolites_file, sep='\t')
    smiles_to_spec_names = {name_to_smiles_dict[idx]:idx for idx in spec_df.index}
    background_molecule_names = [smiles_to_spec_names[s] for s in background_molecules_df['smiles'] \
                                  if s in smiles_to_spec_names.keys() and smiles_to_spec_names[s] in spec_df.index]
    background_specs = spec_df.loc[background_molecule_names,:]
    background_specs = background_specs[~background_specs.index.duplicated(keep='first')]
    return background_specs


def get_rank(array):
    """
    Assign each element in an array its rank in a sorted list
    """
    argsorted = np.flip(np.argsort(array))
    
    ranks = np.zeros(len(array))
    for i, r in enumerate(argsorted):
        ranks[r] = i
    return ranks

def get_ksi (mat):
    """
    mat (np.array): matrix with shape JxK where J is the number of wavelengths and k is 
        the number of spectra
    """
    pinvS = np.linalg.pinv(mat)
    ksi = 1/(np.sqrt(np.sum(mat**2, axis=0) * np.sum( (pinvS**2).transpose(), axis=0)))
    return ksi

# def get_ksi_pairwise (mat,eps = 1e-8):
#     """
#     mat (np.array): matrix with shape Jx2 where J is the number of wavelengths and k is 
#         the number of spectra
#     """
#     # if there is no overlap between spectra
#     if np.sum(mat[0,:] * mat[1,:]) < eps:
#         return 1
#     else:
#         pinvS = np.linalg.pinv(mat)
#         ksi = 1/(np.sqrt(np.sum(mat**2, axis=0) * np.sum( (pinvS**2).transpose(), axis=0)))
#     return ksi

def get_ksi_pairwise (arr1, arr2):
    """
    arr1 (np.array): array with shape 1xJ where J is the number of wavelengths
    arr2 (np.array): array with shape 1xJ where J is the number of wavelengths
    """
    mat = np.stack([arr1, arr2]).transpose()
    pinvS = np.linalg.pinv(mat)
    ksi = 1/(np.sqrt(np.sum(mat**2, axis=0) * np.sum( (pinvS**2).transpose(), axis=0)))
    return ksi[-1]

def get_lboz_w_background(spec, background_specs, normalize=True):    
    """
    spec (np.array):
    background_specs (np.array): 
    """
    
    if normalize:
        spec = spec/np.max(spec)
    mat = np.vstack([background_specs, np.reshape(spec, [1,len(spec)])]).transpose()
#     pinvS = np.linalg.pinv(mat)
    
    ksi = get_ksi (mat)
    return ksi[-1]

def sliding_window_average(arr, window_size):
    """
    Returns a sliding window averaged copy of array
    """
    window_sum = np.cumsum(np.pad(arr, window_size//2, 'edge'))

    window_sum[window_size:] = window_sum[window_size:] - window_sum[:-window_size]
    window_average = window_sum[window_size - 1:] / window_size

    return window_average

def clean_spectra_for_dist_calcs(specs, wl_lb=350, wl_ub=1000, intensity_thresh = 0.1, normalize=True):
    """
    Truncate spectra, optionally normalize, and filter out spectra with 
    absorbance below a threshold in the region of interest.
    """
    c_idxs = [c for c in specs.columns if c > wl_lb and c < wl_ub]
    if normalize:
        norm_filt_tddft = specs.loc[:, c_idxs] 
        gt_thresh = norm_filt_tddft.sum(axis=1)>intensity_thresh
        norm_filt_tddft = norm_filt_tddft.loc[gt_thresh,:] / specs.loc[gt_thresh,:].max(axis=1).values[:, np.newaxis]
        return norm_filt_tddft
    else:
        filt_tddft = specs.loc[:, c_idxs]
        filt_tddft = specs[(specs.sum(axis=1))>intensity_thresh]
        return filt_tddft
    
    
def parse_spectra(specs, intensity_threshold=0.05, derivative_threshold=0.0005, method_split_method='midpoint'):
    
    original_spectra = specs.values/np.max(specs.values, axis=1)[:,np.newaxis]

# parse the spectra by peak -- can be sensitive to parsing parameters
    failed_to_parse = []
    parsed_spectra = []
    for i in range(original_spectra.shape[0]):
        try:
            parsed_spectra.append(parse_spectrum(original_spectra[i,:],intensity_threshold, derivative_threshold, method_split_method, validate=True))
        except:
            print (i)
            print (sum(original_spectra[i,:]))
            print (sum(filt_tddft.iloc[i,:]))
            plt.plot(original_spectra[i,:])
    assert len(failed_to_parse) == 0
    original_to_parsed_idxs_dict = {}
    counter = 0
    for i in range(len(original_spectra)):
        parsed_idxs = []
        for j in parsed_spectra[i]:
            parsed_idxs.append(counter)
            counter += 1
            original_to_parsed_idxs_dict[i] = parsed_idxs
    parsed_spectra = [x for y in parsed_spectra for x in y]
    return parsed_spectra, original_spectra, original_to_parsed_idxs_dict

def parallel_cdist (X, Y, metric, n_jobs=16):
    """
    Compute pairwise distances between two matrices in parallel
    """
    combos = itertools.product(X, Y)
    distances = Parallel(n_jobs=n_jobs)(delayed(metric)(x,y) for x,y in tqdm(combos, total=len(X)*len(Y)))

    return np.reshape(distances, (len(X), len(Y)))
    
def get_uniqueness (specs, dist_metric, parse_spectra_peaks = False,
                    wl_lb=350, wl_ub=1000, normalize=True, topk_frac=0.1,
                   intensity_thresh=0.1, parsing_derivative_thresh = 0.0005, parsing_intensity_thresh=0.05, n_jobs=1, save_dists=None):
    print ("Filtering spectrum dataframe")
    print ("Original dataframe dimensions:", specs.shape)
    specs = clean_spectra_for_dist_calcs(specs, wl_lb=wl_lb, wl_ub=wl_ub, 
                                         intensity_thresh=intensity_thresh, 
                                         normalize=normalize)
    print ("After filtering:", specs.shape)
    
    if parse_spectra_peaks:
        print ('parsing spectral peaks')
        parsed_spectra, original_spectra, original_to_parsed_idxs_dict = parse_spectra(specs, derivative_threshold=parsing_derivative_thresh, intensity_threshold=parsing_intensity_thresh)
        distances = cdist(specs.values, parsed_spectra, metric=dist_metric)
        scores = (specs.index, compute_scores_from_subset(distances, original_spectra, original_to_parsed_idxs_dict, topk_frac))
    else:
        print ('scoring whole spectra')
        if n_jobs == 1:
            distances = cdist(specs.values, specs.values, metric=dist_metric)
        else:
            distances = parallel_cdist(specs.values, specs.values, metric=dist_metric, n_jobs=n_jobs)
        if save_dists is not None:
            np.save(save_dists, distances)
        scores  = (specs.index, compute_scores_from_subset_no_parsing(distances, topk_frac)[0])
    return scores

def interpolate_spectra (spec_df, target_wavelengths):
    interp_df = spec_df.T.copy()
    interp_df['wavelength'] = interp_df.index
    interp_df = interp_df.merge(pd.DataFrame({'wavelength':target_wavelengths}), on='wavelength', how='outer').sort_values(by='wavelength')

    wavelengths = interp_df['wavelength']
    interp_df = interp_df.drop(columns=['wavelength'])
    interp_df.index = wavelengths
    interp_df = interp_df.interpolate(method='index', axis=0)
    interp_df = interp_df.loc[target_wavelengths,:].T
    return interp_df


def get_background_uniqueness (specs, bg_spectra_path, bg_wls_path, pixel_sample_number = 1000,
                    wl_lb=350, wl_ub=1000, normalize=True, intensity_thresh=0.0, intensity_factor=0.1):
    print ("Filtering spectrum dataframe")
    print ("Original dataframe dimensions:", specs.shape)
    #trim and normalize absorbance spectra
    specs = clean_spectra_for_dist_calcs(specs, wl_lb=wl_lb, wl_ub=wl_ub, 
                                         intensity_thresh=intensity_thresh, 
                                         normalize=normalize)
    # read in background spectra
    if isinstance(bg_wls_path,str):
        bg_wls = np.load(bg_wls_path)
    elif isinstance(bg_wls_path, np.ndarray):
        bg_wls = bg_wls_path
    else:
        print ('Format of `bg_wls_path` not currently supported')
        
    # read in background wavelengths
    print ("After filtering:", specs.shape)
    print ("Interpolating to match spectrum and background wavelengths")
    #interpolate absorance spectra to match reflectance
    specs = interpolate_spectra(specs, bg_wls)
    print ("After interpolating:", specs.shape)
    
    
    if isinstance(bg_spectra_path,str):
        bg_specs = np.load(bg_spectra_path)
    elif isinstance(bg_wls_path, np.ndarray):
        bg_specs = bg_spectra_path
    else:
        print ('Format of `bg_spectra_path` not currently supported')
    if pixel_sample_number > bg_specs.shape[0]:
        print ('using all pixels')
        pixel_sample_number = bg_specs.shape[0]
    
    rand_subset_idxs = np.random.choice(range(bg_specs.shape[0]), pixel_sample_number, replace=False)
    bg_specs = bg_specs[rand_subset_idxs, :]
    bg_specs = bg_specs / np.max(bg_specs, axis=1)[:, np.newaxis]
    
    
    scores = get_impact_scores(specs.values, bg_specs, intensity_factor=intensity_factor)
    
    return specs.index, scores

def wasserstein_dist (xs, d1, d2):
    try:
        return stats.wasserstein_distance(xs, xs, d1, d2)
    except ValueError:
        return None

def compare_spectrum_subset (spectrum, subset_to_compare, metric):
    return [metric(spectrum,s) for s in subset_to_compare]

def calc_dist_to_subset (spectra_to_score, subset_to_compare, metric, n_jobs = 16):
    dist_mat = Parallel(n_jobs=n_jobs)(delayed(compare_spectrum_subset)(s, subset_to_compare, metric) for s in tqdm(spectra_to_score))
    return dist_mat