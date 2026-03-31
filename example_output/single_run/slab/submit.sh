#!/bin/bash

#SBATCH -p haswell
#SBATCH -t 30:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=30
#SBATCH -J vasp
#SBATCH --qos=nogpu
#SBATCH --job-name="Slab"

#cd 

source /etc/profile.d/zlmod.sh
module load arch/haswell24v2
module load intel-oneapi-mkl/2024.2.2
module load intel-oneapi-mpi/2021.12.1
module load intel/2025.0.0
module load vasp/5.4.4.pl2

ulimit -s unlimited

mpirun vasp_std > out

exit

