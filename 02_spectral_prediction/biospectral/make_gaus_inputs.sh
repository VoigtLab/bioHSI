#!/bin/bash
OPTIMIZE=true
XTB=true
FUNCTIONAL="b3lyp"
BASIS="6-31g(d)"
MEM=16
CPUS=10
IN_EXT="coord"
GAU_OUTPUT_DIR="/nobackup1c/users/itail/"
OVERWRITE=false

while (( $# )); do
  	case "$1" in
		-i|--input-dir)
				INPUT_DIR=$2
		;;
		-o|--output-dir)
				OUTPUT_DIR=$2
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
		--no-xtb)
				XTB=false
		;;
		--no-opt)
				OPTIMIZE=false 
		;;
		--overwrite)
				OVERWRITE=true
		;;
    -*)   echo "Bad option" >&2; exit 1 ;;
    --)   shift; args+=( "$@" ); set -- ;;
    *)    args+=( "$1" ) ;;
  	esac
  	shift
done


if [ "$XTB" = true ]; then
mkdir "$OUTPUT_DIR/xtb_out"
	for FILE in $(ls $INPUT_DIR/*$IN_EXT); do
		tmp=${FILE##*/}
		NAME=${tmp%.$IN_EXT}

		if [ "$OVERWRITE" = false ] && [ -f  "$OUTPUT_DIR/$NAME.com" ]; then
			echo "File $OUTPUT_DIR/$NAME.com already exists. Skipping."
		else		
			XTBOUT="$OUTPUT_DIR/xtb_out/$NAME"
			input="$OUTPUT_DIR/${NAME}_in.xyz"
                        XTBOUTFILE="$XTBOUT.gau"
			
			if  [ "$OVERWRITE" = false ] && [ -f  "$XTBOUTFILE" ]; then
				echo "File $XTBOUTFILE already exists. Skipping."
			else
				expr $(wc -l < $FILE) - 1 > $input
				echo  >> $input
				tail -n +2 $FILE >> $input

				charge=$(head -n 1 $FILE | awk -F" " '{print $1}')
				mult=$(expr $(head -n 1 $FILE | awk -F" " '{print $2}') - 1)
				xtb $input --namespace $XTBOUT --chrg $charge --uhf $mult --opt vtight > $XTBOUT.output
				
				input="$XTBOUT.xtbopt.xyz"
				head -n 1 $FILE > $XTBOUTFILE
				tail -n +3 $input >> $XTBOUTFILE
			fi	
			COMFILE="$OUTPUT_DIR/$NAME.com"

			bash biospectral/coords_to_com.sh -c $XTBOUTFILE -o $COMFILE -n $NAME -f $FUNCTIONAL -b $BASIS --gau-output-dir $GAU_OUTPUT_DIR --mem $MEM --cpus $CPUS --optimize $OPTIMIZE
		fi
	done
else
	for FILE in $(ls $INPUT_DIR/*$IN_EXT); do
		if [ "$OVERWRITE" = false ] && [ -f  "$OUTPUT_DIR/$NAME.com" ]; then
                        echo "File $OUTPUT_DIR/$NAME.com already exists. Skipping."
                else

			tmp=${FILE##*/}
			NAME=${tmp%.$IN_EXT}
			COMFILE="$OUTPUT_DIR/$NAME.com"
			bash biospectral/coords_to_com.sh -c $FILE -o $COMFILE -n $NAME -f $FUNCTIONAL -b $BASIS --gau-output-dir $GAU_OUTPUT_DIR --mem $MEM --cpus $CPUS --optimize $OPTIMIZE
		fi
	done
fi

