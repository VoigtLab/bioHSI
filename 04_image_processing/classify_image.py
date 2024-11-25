from hsi_detect.utils import *
import argparse
from datetime import date
# Plotting parameters
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def parse_arguments():
  parser = argparse.ArgumentParser(description='Classify images based on reference spectrum.')
  parser.add_argument('--reference-spectrum-path', type=str, required=True, help='Path to the reference spectrum file.')
  parser.add_argument('--image-path', type=str, required=True, help='Path to the image hdr file.')
  parser.add_argument('--save-prefix', type=str, required=True, help='Prefix to save the classified image.')
  parser.add_argument('--dist-threshold', type=float, required=True, help='Distance threshold for classification.')
  parser.add_argument('--smoothing-window-size', type=int, default=11, required=False, help='Window for smoothing pixel spectra.')
  parser.add_argument('--reduced-dims', type=int, default=3, required=False, help='Number of dimensions to keep in PCA step of classication.')
  parser.add_argument('--filter-threshold', type=float, required=True, help='Similarity to reference threshold above which to remove endmembers.')
  return parser.parse_args()

if __name__ == "__main__":
  today = date.today()
  date_str = today.strftime('%d%b%Y')
  print ('Date prefix:', date_str)

  args = parse_arguments()
  reference_spectrum_path = args.reference_spectrum_path
  image_path = args.image_path
  save_prefix = args.save_prefix

  #Load image
  lib = envi.open(image_path)
  centers = lib.bands.centers # wavelength band centers
  unnormalized_img = smooth_img_spectrum(lib.load(), window_size=args.smoothing_window_size)
  img = unnormalized_img / np.nanmax(unnormalized_img, axis=2, keepdims=True)

  # Load reference spectrum 
  # Assuming the reference spectrum is an np array (.npy)
  # With the first row containing the wavelengths
  # And second row containing the intensities of absorbance
  reference_file = np.load(reference_spectrum_path)
  reference_wls = reference_file[0,:]
  reference_spec = reference_file[1,:]

  # interpolate to match the image wavelengths
  reference_spec = np.interp(centers, reference_wls, reference_spec)

  # Detection 
  img_flattened= np.reshape(img, (img.shape[0]*img.shape[1], img.shape[2]))
  img_flattened = img_flattened / np.nanmax(img_flattened, axis=1, keepdims=True)

  em_ls = kmeans_hierarchical_extract_endmembers(img, return_cluster_idxs=False, 
                                                             reference_spec=reference_spec, 
                                                             reduced_dims=args.reduced_dims,
                                                             filter_threshold=args.filter_threshold,
                                                             distance_threshold=args.distance_threshold,
                                                             n_init=5)
  
  mat = np.vstack([-reference_spec, em_ls[0]])
  scored_img = UCLS(img_flattened, mat)
  scored_img = np.reshape(scored_img[:,0], img.shape[:2])

  # threshold score
  scored_img[scored_img<0] = 0

  # Visualize classified image

  plt.figure(dpi=500)
  plt.imshow(scored_img, vmin=0, vmax=0.5, cmap='inferno')
  plt.xticks([])
  plt.yticks([])
  plt.box(False)
  plt.savefig(save_prefix+f'{date_str}_classified.png', dpi=500, transparent=True)
  plt.colorbar()
  plt.savefig(save_prefix+f'{date_str}_classified_w_colorbar.pdf', dpi=500, transparent=True)
  plt.show()