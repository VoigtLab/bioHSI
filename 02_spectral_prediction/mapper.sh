#!/bin/bash
  
echo "Name: " $1
echo "SMILES: " $2

export g16root=/home/gridsan/groups/coley_lab/gaussian  # point to g16 installation folder, modify this line for different server
export PATH=$g16root/g16/:$g16root/gv:$PATH
export GAUSS_EXEDIR=$g16root/g16

module load anaconda/2021a
source activate biospectral
export PYTHONPATH="~/biospectral":$PYTHONPATH

# Run the executable
# The input files must have the name (or identifier) of the molecule in the first column
# and the SMILES in the second column
# all string must be in single quotation marks so that they are read in properly
bash predict_optical_properties.sh -i data/processed/all_accessible_metabolites/ -n "$1" -s "$2" -o data/tddft_outputs/all_accessible/ -m 40 -c 12 -f b3lyp
