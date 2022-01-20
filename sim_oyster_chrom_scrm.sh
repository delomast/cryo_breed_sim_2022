#!/bin/bash
# simulate oyster chromosomes with scrm
# single thread
# output file name $5chr.txt
# takes positional arguments in order:
# 1 number of haplotypes
# 2 Ne, 
# 3 mutation rate per base, 
# 4 expected number of recombinations per chromosome per generation,
# 5 output file prefix
# 6 random seed
# 7 chromosome length

# activate conda environment containing scrm
module load miniconda
source activate /project/oyster_gs_sim/scrm

i=0
# calculate 4*Ne*u*L
t=$(awk "BEGIN {print 4*$2*$3*$7}")
# calculate 4*Ne*r
R=$(awk "BEGIN {print 4*$2*$4}")
scrm $1 1 -t $t -r $R $7 -l 200r -seed $6 -transpose-segsites > "$5"chr.txt
