#!/bin/bash
# run simulations with scrm for
# performance of different SNP panel sizes
# written for execution on Ceres
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=20G # Maximum amount of memory this job will be given
#SBATCH --time=24:00:00 # ask that the job be allowed to run for 
#SBATCH --array=1-2%2 #specify how many jobs in the array and limit number running concurrently (e.g. 1-96%40)
#SBATCH --output=arrayCryo_%a.out # tell it where to store the output console text

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load r/4.1.2

# check for random seeds
if [ ! -f randSeeds.txt ]; then
    echo "randSeeds.txt does not exist."
	exit 1
fi

# get random seed
x=$(cat randSeeds.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo "My random seed is: " $x

# make temp directory
mkdir /90daydata/oyster_gs_sim/cryo/temp"$SLURM_ARRAY_TASK_ID"

# run simulation WITH cryo
# randomSeed iterationNumber TemporaryLocalStorageDirectory useCryo
Rscript scrm_HPC.R $x $SLURM_ARRAY_TASK_ID /90daydata/oyster_gs_sim/cryo/temp"$SLURM_ARRAY_TASK_ID"/ TRUE
# and WITHOUT
Rscript scrm_HPC.R $x $SLURM_ARRAY_TASK_ID /90daydata/oyster_gs_sim/cryo/temp"$SLURM_ARRAY_TASK_ID"/ FALSE

# remove temp directory
rm -r /90daydata/oyster_gs_sim/cryo/temp"$SLURM_ARRAY_TASK_ID"

echo "Done with simulation"
