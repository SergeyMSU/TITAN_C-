#! /bin/bash
#SBATCH --job-name="TITAN_5"
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=0-2:30
#SBATCH --mail-user=sergey.korolkov@cosmos.ru
#SBATCH --mail-type=ALL
#SBATCH --constraint="type_e"
#SBATCH --partition=rocky

# Executable
srun ./my_program