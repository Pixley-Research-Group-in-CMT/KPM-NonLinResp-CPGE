#!/bin/bash

#SBATCH --partition=gpu           # Partition (job queue)

#SBATCH --gpus=1 -c 1             # first 1 GPUs second 1 CPUs

#SBATCH --requeue                 # Return job to the queue if preempted

#SBATCH --job-name=Cg04Ws       # Assign a short name to your job

#SBATCH --nodes=1                 # Number of nodes you require

#SBATCH --ntasks-per-node=1       # Total # of tasks across all nodes

#SBATCH --mem=80000                # Real memory (RAM) required (MB)

#SBATCH --constraint=a100,ib             #a100-80gb,ib# constraint type of node

#SBATCH --time=24:00:00           # Total run time limit in minutes

#SBATCH --array=1-20:1                  # Use job array ID as a variable 

#SBATCH --error=slurm.%N.%j.err   # STDERR output file (optional)
###SBATCH --output=slurm.%N.%j.out  # STDOUT output file

##SBATCH --mail-type=begin           # send an email when the program begins

#SBATCH --mail-type=end,fail           # send an email when the program ends

#SBATCH --mail-user=awu@flatironinstitute.org


mkdir -p $HOME/CPGE/$SLURM_JOB_ID

cd $HOME/CPGE/$SLURM_JOB_ID

###export JULIA_NUM_THREADS=1
srun julia $HOME/CPGE/CPGEfromGammaNCs.jl $SLURM_ARRAY_TASK_ID @> output${SLURM_ARRAY_TASK_ID}.log






