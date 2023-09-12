#!/bin/bash
#SBATCH -s serial
#SBATCH -n 12
#SBATCH -a 1-1
#SBATCH -t 2-00:00:00

./cpptamcmc -S $SLURM_ARRAY_TASK_ID -L $SLURM_ARRAY_TASK_ID
