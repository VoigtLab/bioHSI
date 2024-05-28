#!/bin/bash
TASKNUM=1
while (( $# )); do
    case "$1" in
        -i|--input-dir)
        INPUT_DIR=$2
    ;;
        -o|--output-dir)
        OUTPUT_DIR=$2
    ;;
        -j|--jobs)
        TASKNUM=$2
    ;;
    esac
    shift
done

mkdir $OUTPUT_DIR

FILELIST=$(ls $INPUT_DIR/*com)
FILENUM=$(echo $FILELIST |  wc -w)

BATCHSIZE=$(expr $FILENUM / $TASKNUM + 1)

for (( STARTLINE=0; STARTLINE<$FILENUM ; STARTLINE+=$BATCHSIZE )); do
    INPUTS=$(ls $INPUT_DIR/*com | tail -n +${STARTLINE} | head -n $BATCHSIZE)
    NAME="$STARTLINE"_inputs
    echo "$OUTPUT_DIR/$NAME.err"
    sbatch -J $NAME -e $OUTPUT_DIR/25May2022_$NAME.err -o $OUTPUT_DIR/25May2022_$NAME.out biospectral/run_script.sh -i "$INPUTS" -o $OUTPUT_DIR
done
