# biospectral

Predict absorbance spectra for a set of SMILES

### 1. Generate Gaussian 16 input files 
`db_to_gaussian_input.sh`\
Parameters:\
--smiles-file: `string` CSV file containing input SMILES\
--name-col: `string` name of column in smiles-file containing name (unique IDs) of molcules used to label output files\
--smiles-col: `string` name of column in smiles-file containing SMILES\
-d, --dir: `string` directory to store generated 3D structure files and Gaussian input files\
-f, --functional: `string` of Gaussian 16 functional to be used in TD-DFT calculation\
-b, --basis: `string` of Gaussian 16 basis set to be used in TD-DFT calculation\
--mem: `int` memory to be allocated for TD-DFT jobs\
--cpus: `int` cpus to be allocated for TD-DFT jobs
 
### 2. Predict absorbance spectra using TD-DFT with Gaussian 16
`tddft_batched.sh`\
Parameters:\
--input-dir: `string` directory containing Gaussian 16 input files\
--output-dir: `string` directory for saving output files\
--jobs: `int` number of separate Slurm jobs to split the work into

### 2.5 Retrieve singlet excitation values from Gaussian 
run    `grep Singlet *log` in the output directory of Gaussian predictions\
TODO: incorporate this with the rest of the workflow

### 3. Generate UV-Vis absorbance spectra from TD-DFT Singlets
`get_spectra.py`\
Parameters:\
--singlets: `string` file name with excitation information from TD-DFT predictions (generated in step 2.5)\
--out-file: `string` name for CSV file in which spectra are saved\
--lower-bound: `int` lower bound (in nm) for the wavelengths for spectrum generation\
--upper-bound: `int` upper bound (in nm) for the wavelengths for spectrum generation\
-s: `float` parameter determining gaussian FWHM during convolution
