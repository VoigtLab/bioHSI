#!bin/bash

#The variable TMPDIR may need to be set

# Runs the pipeline from SMILES through the TD-DFT simulation for one molecule.

while getopts i:n:s:o:m:c:f: option
	do
		case "${option}"
			in
					i) INPUT_DIR=${OPTARG};;
					n) NAME=${OPTARG};;
          s) SMILES=${OPTARG};;
					o) OUTPUT_DIR=${OPTARG};;
					m) MEM=${OPTARG};;
					c) CPUS=${OPTARG};;
					f) FUNCTIONAL=${OPTARG};;
		esac
done


echo "Getting 3D structure from SMILES..."
python biospectral/make_structure_single.py --smiles $SMILES --name $NAME --out-dir $INPUT_DIR
echo "...done getting 3D structure from SMILES"

echo "Refining 3D structure and making Gaussian input file..."
bash biospectral/make_gaus_inputs_single.sh --name $NAME -i $INPUT_DIR -o $INPUT_DIR --gau-output-dir $TMPDIR --mem $MEM --cpus $CPUS
echo "...done refining 3D structure and making Gaussian input file"

echo "Redefining Gaussian input file for supercloud..."
bash biospectral/make_diff_functional_inputs_single.sh --name $NAME -i $INPUT_DIR -o $INPUT_DIR --gau-output-dir $TMPDIR --mem $MEM --cpus $CPUS -f $FUNCTIONAL
echo "...done redefining Gaussian input file for supercloud"

echo "Running TD-DFT simulation..."
bash biospectral/run_script_single.sh -i $INPUT_DIR -n $NAME"_$FUNCTIONAL" -o $OUTPUT_DIR
echo "...done running TD-DFT simulation"