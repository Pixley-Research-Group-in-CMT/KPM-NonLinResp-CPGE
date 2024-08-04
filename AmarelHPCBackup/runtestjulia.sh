#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)

#SBATCH --requeue                 # Return job to the queue if preempted

#SBATCH --job-name=test       # Assign a short name to your job

#SBATCH --nodes=1                 # Number of nodes you require

#SBATCH --ntasks=1                # Total # of tasks across all nodes

#SBATCH --cpus-per-task=10         # cpu-Cores per task (>1 if multithread tasks)

#SBATCH --mem=40000                # Real memory (RAM) required (MB)

#SBATCH --time=24:00:00           # Total run time limit (HH:MM:SS)

#SBATCH --output=slurm.%N.%j.out  # STDOUT output file

#SBATCH --error=slurm.%N.%j.err   # STDERR output file (optional)

#SBATCH --mail-type=begin         # send an email when the program begins

#SBATCH --mail-type=end           # send an email when the program ends

#SBATCH --mail-user=aw666@scarletmail.rutgers.edu

module purge
module load intel/19.0.3

mkdir -p $HOME/CPGE/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID

cd $HOME/CPGE/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID

export JULIA_NUM_THREADS=2
srun --mpi=pmi2 julia -p 1 $HOME/CPGE/rarestates.jl $SLURM_ARRAY_TASK_ID @> output${SLURM_ARRAY_TASK_ID}.log

sleep 3

sacct --format NTasks,MaxRSS,MaxVMSize,AveRSS,AveVMSize,AveCPU -j $SLURM_JOBID

sleep 2
         

