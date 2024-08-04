#!/bin/bash

#SBATCH --partition=gpu          # Partition (job queue)

#SBATCH --requeue                 # Return job to the queue if preempted

#SBATCH --job-name=Gyxz       # Assign a short name to your job

#SBATCH --nodes=1                 # Number of nodes you require

#SBATCH --ntasks=1                # Total # of tasks across all nodes

#SBATCH --cpus-per-task=1         # cpu-Cores per task (>1 if multithread tasks)

#SBATCH --gres=gpu:1             # number of gpus per node

#SBATCH --constraint="volta|pascal|titan" # add constraint to the GPU nodeampere

#SBATCH --mem=64000                # Real memory (RAM) required (MB)

#SBATCH --time=24:00:00           # Total run time limit (HH:MM:SS)

#SBATCH --array=1-100:1           # Use job array ID as a variable

#SBATCH --output=slurm.%N.%j.out  # STDOUT output file

#SBATCH --error=slurm.%N.%j.err   # STDERR output file (optional)

#SBATCH --mail-type=begin           # send an email when the program begins

#SBATCH --mail-type=end,fail           # send an email when the program ends

#SBATCH --mail-user=aw666@scarletmail.rutgers.edu

mkdir -p $HOME/CPGE/$SLURM_JOB_ID
cd $HOME/CPGE/$SLURM_JOB_ID

###export JULIA_NUM_THREADS=1
srun julia $HOME/CPGE/collectGammaPartial.jl $SLURM_ARRAY_TASK_ID @> output${SLURM_ARRAY_TASK_ID}.log
#srun julia $HOME/CPGE/collectKPMmu.jl $SLURM_ARRAY_TASK_ID @> output${SLURM_ARRAY_TASK_ID}.log

##sleep 3

##sacct --format NTasks,MaxRSS,MaxVMSize,AveRSS,AveVMSize,AveCPU -j $SLURM_JOBID

##sleep 2


