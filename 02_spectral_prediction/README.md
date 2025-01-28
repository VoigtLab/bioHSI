# biospectral

Predict absorbance spectra for a set of SMILES using Gaussian 16. Running the TD-DFT predictions requires a Gaussian 16 license and installation. 

To run the scripts that generate the Gaussian 16 inputs and call the TD-DFT prediction code, scripts are provided to submit jobs via various job scheduler interfaces (e.g., SLURM, and a system unique to MIT Lincoln Lab's Supercloud cluster). The underlying pipeline is the same. 

Note that the TD-DFT simulations take a long time to run and can take months to run for 1000s of molecules on multiple CPUs. 


## Step 1A. Run through SLURM in batched mode
### A1. Generate Gaussian 16 input files 
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
 
### A2. Predict absorbance spectra using TD-DFT with Gaussian 16
`tddft_batched.sh`\
Parameters:\
--input-dir: `string` directory containing Gaussian 16 input files\
--output-dir: `string` directory for saving output files\
--jobs: `int` number of separate Slurm jobs to split the work into

## Step 1B. Run single jobs 
(can be run in batch on MIT Lincoln Lab Supercloud cluster using `LLMapReduce --mapper mapper.sh --input <path_to_raw_data> --output output/ --np <cpu usage specification>`)

`predict_optical_properties.sh`\
Parameters:\
-i: `string` Path to directory for saving generated 3D molecular coordinates files\
-n: `string` Name of input molecule (used as label for saving files)\
-s: `string` SMILES string of input molecule\
-o: `string` Path to directory for saving TD-DFT prediction outputs\
-m: `int` Memory to allocate for Gaussian 16 TD-DFT job\
-c: `int` CPU number to allocate for Gaussian 16 TD-DFT job\
-f: `string` Name of TD-DFT functional to use for Gaussian 16 TD-DFT job.\

### Step 2.5 Retrieve singlet excitation values from Gaussian 
run    `grep Singlet *log` in the output directory of Gaussian predictions\
TODO: incorporate this with the rest of the workflow

### Step 3. Generate UV-Vis absorbance spectra from TD-DFT Singlets
`get_spectra.py`\
Parameters:\
--singlets: `string` file name with excitation information from TD-DFT predictions (generated in step 2.5)\
--out-file: `string` name for CSV file in which spectra are saved\
--lower-bound: `int` lower bound (in nm) for the wavelengths for spectrum generation\
--upper-bound: `int` upper bound (in nm) for the wavelengths for spectrum generation\
-s: `float` parameter determining gaussian FWHM during convolution\
--fixed-width: if flag is set, all peaks have equal width with FWHM determined b s in nm.\