#!bin/bash

while getopts i:n:s:o:m:c: option
	do
		case "${option}"
			in
					i) INPUT_DIR=${OPTARG};;
					n) NAME=${OPTARG};;
          				s) SMILES=${OPTARG};;
					o) OUTPUT_DIR=${OPTARG};;
					m) MEM=${OPTARG};;
					c) CPUS=${OPTARG};;
		esac
done


echo "Getting 3D structure from SMILES..."
python biospectral/make_structure_single.py --smiles $SMILES --name $NAME --out-dir $INPUT_DIR
echo "...done getting 3D structure from SMILES"

echo "Refining 3D structure and making Gaussian input file..."
bash biospectral/make_gaus_inputs_single.sh --name $NAME -i $INPUT_DIR -o $INPUT_DIR --gau-output-dir $TMPDIR --mem 40 --cpus 10
echo "...done refining 3D structure and making Gaussian input file"

echo "Redefining Gaussian input file for supercloud..."
bash biospectral/make_diff_functional_inputs_single.sh --name $NAME -i $INPUT_DIR -o $INPUT_DIR --gau-output-dir "/home/gridsan/itail/scratch/gaussian/" --mem $MEM --cpus $CPUS -f "b3lyp"
echo "...done redefining Gaussian input file for supercloud"

echo "Running TD-DFT simulation..."
bash biospectral/run_script_single.sh -i $INPUT_DIR -n $NAME"_b3lyp" -o $OUTPUT_DIR
echo "...done running TD-DFT simulation"


