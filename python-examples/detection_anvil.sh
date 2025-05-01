#!/bin/bash
#SBATCH --job-name=detect
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH -A cis220161
#SBATCH -p wholenode
#SBATCH -o /anvil/projects/x-cis220161/datasets/aps/apr17-metrics/sz3detect-%j.o      # Name of stdout output file
#SBATCH -e /anvil/projects/x-cis220161/datasets/aps/apr17-metrics/sz3detect-%j.e      # Name of stderr error file

module load anaconda
conda activate rare_event
cd /home/x-yliu4/aps/SZ3/python-examples
python -u detection_benchmark.py