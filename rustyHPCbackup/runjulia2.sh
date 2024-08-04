#!/bin/bash
#SBATCH --partition=ccq               # Partition (job queue) gen/ccq
###SBATCH --gpus=1 -c 1                  # first 1 GPUs second 1 CPUs
#SBATCH --requeue                      # Return job to the queue if preempted
#SBATCH --job-name=H0gam01         # Assign a short name for the job
#SBATCH --nodes=1                      # Number of nodes you require
#SBATCH --ntasks-per-node=1            # Total number of tasks across all nodes
#SBATCH --cpus-per-task=20             # cpu-cores per task (>1 if multithread tasks)
####SBATCH --constraint=rome            # Constraint the type of node broadwell/skylake/rome
#SBATCH --mem=1000000                    # Real memory (RAM) required (MB)
#SBATCH --time=24:00:00                # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out       # STDOUT output file
#SBATCH --error=slurm.%N.%j.err        # STDOUT output file (optional)
##SBATCH --mail-type=begin              # send an email when the program begins
#SBATCH --mail-type=end,fail                # send an email when it ends
#SBATCH --mail-user=awu@flatironinstitute.org

mkdir -p $HOME/CPGE/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID

cd $HOME/CPGE/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
export JULIA_NUM_THREADS=1
srun julia $HOME/CPGE/storeHam0.jl $SLURM_ARRAY_TASK_ID @> output${SLURM_ARRAY_TASK_ID}.log

