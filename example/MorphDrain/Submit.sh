#!/bin/bash
#SBATCH --job-name=SP
#SBATCH -N 1
#SBATCH -p debug_queue
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00


module load openmpi_3.0.0/gcc_6.3.0
module load gcc/6.3.0

#date

mpirun --bind-to core:overload-allowed $HOME/LBPM_Master_Chris_Build/tests/RayTrace2 input.db

date

mpirun --bind-to core:overload-allowed $HOME/LBPM_Master_Chris_Build/tests/lbpm_random_force_simulator input.db

date
