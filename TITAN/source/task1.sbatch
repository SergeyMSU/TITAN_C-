#! /bin/bash
#SBATCH --job-name="Hello"
#SBATCH --gpus=0
#SBATCH --nodes=1                      # ���������� ������������ �����
#SBATCH --ntasks=1                     # ���������� MPI ���������
#SBATCH --cpus-per-task=44             # ���������� ���� CPU, ���������� ��� ������ ��������
#SBATCH --cores-per-socket=22          # ���������� ���� �� ����� ������ ����
#SBATCH --time=0-0:20
#SBATCH --mail-user=sergey.korolkov@cosmos.ru
#SBATCH --mail-type=ALL
#SBATCH --constraint="type_b"

# Executable
srun ./hello
