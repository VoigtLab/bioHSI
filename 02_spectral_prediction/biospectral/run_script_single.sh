while getopts i:n:o: option
	do
		case "${option}"
			in
					i) INPUT_DIR=${OPTARG};;
					n) NAME=${OPTARG};;
					o) OUTPUT_DIR=${OPTARG};;
		esac
done

date

INPUT=$INPUT_DIR/$NAME.com
OUTPUT=$OUTPUT_DIR/$NAME.log
echo $INPUT
echo $OUTPUT
echo $(which g16)

if [ -f  "$OUTPUT" ] && [ $(grep Singlet $OUTPUT | wc -l) -gt 0 ]; then
		echo "Absorbance for $NAME already predicted. Skipping."
else
	g16 < $INPUT > $OUTPUT
fi



echo "Job done! Calculation for $NAME completed."
date
