#!/bin/bash
##SBATCH -p preempt_std
##SBATCH -t 0:20:00
#SBATCH -t 1-23:59:59
#SBATCH -p serial
#SBATCH -n 10
#SBATCH -C avx2
#SBATCH -a 1-44

module purge
module load gcc/9.2.0
module load cmake/3.18.4
module load gsl/2.5
#module load eigen/3.3.8
module load boost/1.74.0
module load gnuplot/5.4.3

#export EIGEN3_INCLUDE_DIR=$EIGEN_INC/eigen3
export EIGEN3_INCLUDE_DIR=/home/ob19/Programs/eigen-3.4.0
export CXX=g++
export OMP_NUM_THREADS=10

echo "Processing star position=" $SLURM_ARRAY_TASK_ID

./cpptamcmc execute 1 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID
