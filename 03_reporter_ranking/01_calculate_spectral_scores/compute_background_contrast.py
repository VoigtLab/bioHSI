import pandas as pd
from spectral import *
from scipy.spatial.distance import pdist, cdist, squareform, cosine, jensenshannon
import numpy as np
from spectranalysis.utils import *
from spectranalysis.uniqueness_utils import *
from datetime import date
from joblib import Parallel, delayed
import argparse

# sample submission 
# python compute_background_impact.py --input_spectra background_spectra/pixels_from_arm.npy\ 
# --input_wavelengths background_spectra/wavelengths_from_arm.npy --output_file background_spectra/arm_uniqueness.npy\ 
# --method wasserstein --n_jobs 16

def parse_args():
    parser = argparse.ArgumentParser(description='Compute background impact scores for spectra')
    parser.add_argument('--input_background', type=str, help='Path to input spectra file')
    parser.add_argument('--output_file', type=str, help='Path to output file')
    parser.add_argument('--n_jobs', type=int, default=16, help='Number of jobs to run in parallel')
    return parser.parse_args()

def get_distance_score_ (spectrum, bg_spectra, metric, aggregate=np.nanmean):
    """
    score one spectrum against a set of background spectra using the designated distance metric
    """
    dists = [metric(spectrum, bg) for bg in bg_spectra]
    if any (dists):
        return aggregate(dists)
    else:
        return None

def get_distance_scores(spectra_to_score, bg_spectra, metric, n_jobs=16, aggregate=np.nanmean):
    distance_scores = Parallel(n_jobs=n_jobs)(delayed(get_distance_score_)(s, bg_spectra, metric, aggregate) for s in tqdm(spectra_to_score))
    return distance_scores

def get_diff_pixels_from_image(pixels, n_samples=10000, normalize=True, seed=55):
    np.random.seed(seed)
    idx_1 = np.random.choice(range(pixels.shape[0]), n_samples, replace=True)
    idx_2 = np.random.choice(range(pixels.shape[0]), n_samples, replace=True)
    if normalize:
        diff_pixels = pixels[idx_1,:]/np.max(pixels[idx_1,:], axis=1, keepdims=True) - pixels[idx_2,:]/np.max(pixels[idx_2,:], axis=1, keepdims=True)
    else:
        diff_pixels = pixels[idx_1,:] - pixels[idx_2,:]
    return diff_pixels

def _calc_contrast_to_diff (spec, diff_pixels):
    diff_vals = np.zeros(len(diff_pixels))
    for i,p in enumerate(diff_pixels):
        diff_vals[i] = 1-cosine(p, spec)
    score = np.mean(np.abs(diff_vals))
    return score

def calc_contrast_to_diffs(specs, diff_pixels, n_jobs=16):
    scores = Parallel(n_jobs=n_jobs)(delayed(_calc_contrast_to_diff)(s, diff_pixels) for s in tqdm(specs))
    return scores

if __name__=='__main__':
    args = parse_args()
    all_tddft = pd.read_csv('input_data/16Jul2023_all_done_spectra.csv', index_col=0).transpose()
    all_tddft.index = [m.split('/')[-1].replace('_b3lyp','') for m in all_tddft.index]

    # Set bounds for wavelengths to compare (in nm)
    WL_LOWER_BOUND = 0
    WL_UPPER_BOUND = 5000
    NUM_SPECTRA = 1000
    INTENSITY_FACTOR = 0.1

    # Set minimum intensity threhold for predicted spectra
    # Spectra with a sum of values below this value are filtered out
    INT_THRESHOLD = 0.0
    all_tddft.columns = all_tddft.columns.astype(float)
    all_tddft = all_tddft.loc[~all_tddft.index.duplicated(), :]

    c_idxs = [c for c in all_tddft.columns if c > WL_LOWER_BOUND and c < WL_UPPER_BOUND]


    # normalize?
    norm_filt_tddft = all_tddft.loc[:, c_idxs] 
    gt_thresh = norm_filt_tddft.sum(axis=1)>INT_THRESHOLD
    norm_filt_tddft = norm_filt_tddft.loc[gt_thresh,:] / all_tddft.loc[gt_thresh,:].max(axis=1).values[:, np.newaxis]
    
    # norm_filt_tddft = norm_filt_tddft.loc[gt_thresh,:] / np.trapz(all_tddft.loc[gt_thresh,:].values, x=c_idxs, axis=1)[:, np.newaxis]
    lib = envi.open(args.input_background)
    img = np.array(lib.load())

    pixels  = []
    pixel_idxs = []
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if any(img[i, j, :]):
                pixels.append(img[i, j, :])
                pixel_idxs.append((i,j))

    spectra = np.stack(pixels)
    wls = lib.bands.centers

    interp_df = norm_filt_tddft.T.copy()
    interp_df['wavelength'] = interp_df.index
    interp_df = interp_df.merge(pd.DataFrame({'wavelength':wls}), on='wavelength', how='outer').sort_values(by='wavelength')

    wavelengths = interp_df['wavelength']
    interp_df = interp_df.drop(columns=['wavelength'])
    interp_df.index = wavelengths
    interp_df = interp_df.interpolate(method='index', axis=0)

    interp_df = interp_df.loc[wls,:].T

    

    # distance to differences
    predicted_spectra = interp_df.values
    difference_pixels = get_diff_pixels_from_image(spectra, n_samples=10000, normalize=True)
    np.save(args.output_file+'_difference_pixels.npy', difference_pixels)
    bg_impact = calc_contrast_to_diffs(predicted_spectra, difference_pixels, n_jobs=args.n_jobs)




    output = np.array([interp_df.index, bg_impact])
    np.save(args.output_file, output)
