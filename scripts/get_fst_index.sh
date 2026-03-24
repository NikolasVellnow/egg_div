#!/bin/bash -l


# type in file name 1 without "saf.idx" ending (to provide saf files)
NAME1=$1

# type in file name 2 without "saf.idx" ending (to provide saf files)
NAME2=$2

# number of threads used in angsd
NUM_THREADS=$3

T0=$(date +%T)
echo "Start data processing:"
echo $T0


# Properly initialize conda for non-interactive shells
source ~/miniconda3/etc/profile.d/conda.sh

conda activate angsd

realSFS \
fst \
index \
"${NAME1}".saf.idx \
"${NAME2}".saf.idx \
-sfs "${NAME1}"_"${NAME2}".ml \
-fold 1 \
-P $NUM_THREADS \
-fstout "${NAME1}"_"${NAME2}"


conda deactivate

T1=$(date +%T)
echo "Finished data processing:"
echo $T1

echo "done"

