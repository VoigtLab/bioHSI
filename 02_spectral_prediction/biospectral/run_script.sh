#!/bin/bash 
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --cpus-per-task=5
#SBATCH --time=14-00:00:00
#SBATCH --partition=sched_mit_ccoley
#SBATCH --nodelist=node1238
#SBATCH --mem=18000 

while getopts i:o: option
do
case "${option}"
in
		i) INPUTS=${OPTARG};;
		o) OUTPUT_DIR=${OPTARG};;
esac
done

module add gaussian/16.c01_avx2
date
for INPUT in $(echo $INPUTS | tr " " "\n"); do 
		echo $INPUT
		tmp=${INPUT##*/}
		NAME=${tmp%.com}
		OUTPUT=$OUTPUT_DIR/$NAME.log

		g16 < $INPUT > $OUTPUT

		echo "Job done! Calculation for $NAME completed."
		date

done
