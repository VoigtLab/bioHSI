#!/bin/bash
OPTIMIZE=true
XTB=true
BASIS="6-31g(d)"
MEM=16
CPUS=10
IN_EXT="gau"
GAU_OUTPUT_DIR="/nobackup1c/users/itail/"
while (( $# )); do
  case "$1" in
    -i|--input-dir)
                INPUT_DIR=$2
                ;;
    -o|--output-dir)
                OUTPUT_DIR=$2
                ;;
    --optimize)
                OPTIMIZE=$2
                ;;
    -f|--functionals)
                FUNCTIONALS=$2
                ;;
    -b|--basis)
                BASIS=$2
                ;;
    --gau-output-dir)
                GAU_OUTPUT_DIR=$2
                ;;
    --extension)
				IN_EXT=$2
				;;
    --mem)
                MEM=$2
                ;;
    --cpus)
                CPUS=$2
                ;;
    --name)
                NAME=$2
                ;;
    -*)   echo "Bad option $1" >&2; exit 1 ;;
    --)   shift; args+=( "$@" ); set -- ;;
    *)    args+=( "$1" ) ;;
  esac
  shift
done


for functional in $(echo "$FUNCTIONALS" | tr " " "\n"); do
		FILE="$INPUT_DIR/xtb_out/$NAME.$IN_EXT"
    OUTFILE="$OUTPUT_DIR/"$NAME"_"$functional".com"
		bash biospectral/coords_to_com.sh -c $FILE -o $OUTFILE -n $NAME -f $functional -b $BASIS --gau-output-dir $GAU_OUTPUT_DIR --mem $MEM --cpus $CPUS --optimize $OPTIMIZE
done

