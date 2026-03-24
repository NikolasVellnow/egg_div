#!/bin/bash -l

# type in path to text file with list of samples
SAMPLE_LIST=$1

# type in output file name
OUT=$2

# number of threads used in angsd
NUM_THREADS=$3

# path to reference genome
PATH_REF=/home/mnikvell/Desktop/work/data/genomes/refseq/vertebrate_other/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna


T0=$(date +%T)
echo "Start data processing:"
echo $T0

# Properly initialize conda for non-interactive shells
source ~/miniconda3/etc/profile.d/conda.sh

conda activate angsd

angsd \
-b $SAMPLE_LIST \
-r NC_031775.1 \
-doSaf 1 \
-out ${OUT} \
-anc $PATH_REF \
-GL 2 \
-minMapQ 30 \
-minQ 20 \
-P $NUM_THREADS


conda deactivate

T1=$(date +%T)
echo "Finished data processing:"
echo $T1

echo "done"

