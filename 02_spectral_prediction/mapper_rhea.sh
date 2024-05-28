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
bash predict_optical_properties.sh -i data/processed/rhea_metabolites/ -n "$1" -s "$2" -o data/tddft_outputs/rhea_metabolites/ -m 40 -c 12
