#! /bin/bash
#SBATCH --job-name="TITAN_5"
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=0-0:10
#SBATCH --mail-user=sergey.korolkov@cosmos.ru
#SBATCH --mail-type=ALL
#SBATCH --constraint="type_e"
#SBATCH --partition=rocky

# Executable
module load CUDA
module load nvidia_sdk/nvhpc 
module load gnu7   # или любую другую доступную версию GCC
export GCC_HOME=$(dirname $(dirname $(which gcc)))
#srun make clean
srun make