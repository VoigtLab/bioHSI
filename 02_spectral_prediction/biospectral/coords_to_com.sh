#!/bin/bash
while (( $# )); do
  case "$1" in
    -c|--coordinates)
                COORD=$2
                ;;
    -o|--output-file)
                OUTPUT_FILE=$2
                ;;
    -n|--name)
	    	    NAME=$2
		        ;;
    --optimize)
                OPTIMIZE=$2
                ;;
    -f|--functional)
                FUNCTIONAL=$2
                ;;
    -b|--basis)
                BASIS=$2
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
    -*)   echo "Bad option" >&2; exit 1 ;;
    --)   shift; args+=( "$@" ); set -- ;;
    *)    args+=( "$1" ) ;;
  esac
  shift
done

if [ $OPTIMIZE = true ]; then
    echo '%chk=opt_'$NAME > $OUTPUT_FILE
    echo '%RWF='"$GAU_OUTPUT_DIR""$NAME" >> $OUTPUT_FILE
    echo '%Mem='"$MEM"'GB' >> $OUTPUT_FILE
    echo '%NProcShared='$CPUS >> $OUTPUT_FILE
    echo '#'T opt $FUNCTIONAL/$BASIS >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE
    echo "optimize groundstate of $NAME" >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE
    cat $COORD >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE
    echo "--Link1--" >> $OUTPUT_FILE
    echo '%chk=opt_'$NAME >> $OUTPUT_FILE
    echo '%RWF='"$GAU_OUTPUT_DIR""$NAME" >> $OUTPUT_FILE >> $OUTPUT_FILE
    echo '%NoSave' >> $OUTPUT_FILE
    echo '%Mem='$MEM'GB' >> $OUTPUT_FILE
    echo '%NProcShared='$CPUS >> $OUTPUT_FILE
    echo '#T td=(singlets, nstates=5) scrf=(solvent=water)' $FUNCTIONAL/$BASIS 'geom=check guess=read' >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE
    echo "Calculate excited states using TD-DFT for $NAME" >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE
    head -n 1 $COORD  >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE

else 
    echo '%chk=opt_'$NAME >> $OUTPUT_FILE
    echo '%RWF='"$GAU_OUTPUT_DIR""$NAME" >> $OUTPUT_FILE >> $OUTPUT_FILE
    echo '%NoSave' >> $OUTPUT_FILE
    echo '%Mem='$MEM'GB' >> $OUTPUT_FILE
    echo '%NProcShared='$CPUS >> $OUTPUT_FILE
    echo '#T td=(singlets, nstates=5) scrf=(solvent=water)' $FUNCTIONAL/$BASIS 'geom=check guess=read' >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE
    echo "Calculate excited states using TD-DFT for $NAME" >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE
    cat $COORD >> $OUTPUT_FILE
    echo   >> $OUTPUT_FILE

fi


