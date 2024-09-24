import pandas as pd
import numpy as np
from utils import *
from uniqueness_utils import *
from datetime import date
from joblib import Parallel, delayed
from scipy import stats
import os
import argparse

today = date.today()
date_str = today.strftime('%d%b%Y')
print ('Date prefix:', date_str)

def parse_args():
    parser = argparse.ArgumentParser(description='Compute background impact scores for spectra')
    parser.add_argument('--input_file', type=str, help='Path to input file', default = 'input_data/16Jul2023_all_done_spectra.csv')
    parser.add_argument('--output_file', type=str, help='Path to output file')
    parser.add_argument('--n_jobs', type=int, default=16, help='Number of jobs to run in parallel')
    parser.add_argument('--lb', type=int, default=0, help='Spectral lower bound')
    parser.add_argument('--ub', type=int, default=1000, help='Spectral upper bound')
    return parser.parse_args()

if __name__=='__main__':
    args = parse_args()
    input_file = args.input_file
    output_file = args.output_file
    all_tddft = pd.read_csv(input_file, index_col=0).transpose()
    all_tddft.index = [m.split('/')[-1].replace('_b3lyp','').replace('_wb97xd','') for m in all_tddft.index]

    # Set bounds for wavelengths to compare (in nm)
    WL_LOWER_BOUND = args.lb
    WL_UPPER_BOUND = args.ub

    # Set minimum intensity threhold
    # Spectra with a sum of values below this value are filtered out
    INT_THRESHOLD = 0.0 # need to match to default values from clean_spectra_for_dist_calcs
    all_tddft.columns = all_tddft.columns.astype(float)
    all_tddft = all_tddft.loc[~all_tddft.index.duplicated(), :]
    c_idxs = [c for c in all_tddft.columns if c > WL_LOWER_BOUND and c < WL_UPPER_BOUND]
    all_tddft = all_tddft.loc[:, c_idxs] 
    gt_thresh = all_tddft.sum(axis=1)>INT_THRESHOLD
    all_tddft = all_tddft.loc[gt_thresh,:]
    
    path_name, ext = os.path.splitext(output_file)
    dists_path = f'{path_name}_dist{ext}'
    kwargs = {'wl_lb':WL_LOWER_BOUND, 'wl_ub':WL_UPPER_BOUND}
    # distance_function  = lambda x,y : stats.wasserstein_distance (all_tddft.columns.tolist(), all_tddft.columns.tolist(), x,y)
    all_cdfs = (all_tddft/np.trapz(all_tddft.values, axis=0)).cumsum(axis=0)
    x = np.array(all_tddft.columns.tolist())
    distance_function = lambda x, y: np.trapz(np.abs(x-y), x = x)    
    wasserstein_dist = get_uniqueness(all_cdfs, distance_function, parse_spectra_peaks=False, 
                                      topk_frac = 1, n_jobs=args.n_jobs, 
                                      save_dists=dists_path,intensity_thresh=INT_THRESHOLD, **kwargs)

    np.save(output_file, wasserstein_dist)
