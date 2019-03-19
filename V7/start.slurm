#!/bin/sh
#SBATCH --job-name=V7_OMP
#SBATCH --time=12:00:00
#SBATCH --partition=normal
## Pour reserver 1 noeud
#SBATCH --nodes=1

export OMP_NUM_THREADS=8
time ./a.out
