#!/bin/bash
FUNCTIONAL="b3lyp"
BASIS="6-31g(d)"
MEM=16
CPUS=10
IN_EXT="coord"
GAU_OUTPUT_DIR="/nobackup1c/users/itail/"
NAMECOL=name
SMILESCOL=smiles

while (( $# )); do
  	case "$1" in
		-d|--dir)
				INPUT_DIR=$2
		;;
		-f|--functional)
				FUNCTIONAL=$2
		;;
		-b|--basis)
				BASIS=$2
		;;
		--input-extension)
				IN_EXT=$2
		;;
		--gau-output-dir)
				GAU_OUTPUT_DIR=$2
		;;
		--mem)
				MEM=$2
		;;
		--cpus)
				CPUS=$2
		;;
		--smiles-col)
				SMILESCOL=$2
		;;
		--name-col)
				NAMECOL=$2
		;;
		--smiles-file)
				MOLFILE=$2
		;;
		-*)   echo "Bad option" >&2; exit 1 ;;
		--)   shift; args+=( "$@" ); set -- ;;
		*)    args+=( "$1" ) ;;
  	esac
  	shift
done

module add anaconda3/2020.11
source activate biospectral

export PYTHONPATH="~/biospectral":$PYTHONPATH

echo "Running make_structure.py"
~/.conda/envs/biospectral/bin/python biospectral/make_structure.py --mol-file $MOLFILE --smiles-col $SMILESCOL --name-col $NAMECOL --out-dir $INPUT_DIR

echo "Running make_gaus_inputs.sh"
bash biospectral/make_gaus_inputs.sh -i $INPUT_DIR -o $INPUT_DIR -f $FUNCTIONAL -b $BASIS --input-extension $IN_EXT --gau-output-dir $GAU_OUTPUT_DIR --mem $MEM --cpus $CPUS 
