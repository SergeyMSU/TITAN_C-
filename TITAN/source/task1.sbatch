#! /bin/bash
#SBATCH --job-name="Golikov_task"
#SBATCH --gpus=0
#SBATCH --nodes=50
#SBATCH --gpus-per-node=0
#SBATCH --ntasks=1
#SBATCH --time=0-0:20
#SBATCH --mail-user=sergey.korolkov@cosmos.ru
#SBATCH --mail-type=ALL
#SBATCH --constraint="type_e"

# Executable
srun ./hello
