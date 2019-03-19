#!/bin/sh
#SBATCH --job-name=MERCIER_V3
#SBATCH --time=00:50:00
#SBATCH --partition=short
## Pour reserver 1 noeud
#SBATCH --nodes=1

export OMP_NUM_THREADS=4
time ./a.out
