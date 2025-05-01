#!/bin/bash
#SBATCH --job-name=rare_c
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=01:00:00
#SBATCH -A cis220161
#SBATCH -p wholenode
#SBATCH -o /anvil/projects/x-cis220161/datasets/aps/apr17-metrics/rare_compress_logscale-%j.o      # Name of stdout output file
#SBATCH -e /anvil/projects/x-cis220161/datasets/aps/apr17-metrics/rare_compress_logscale-%j.e      # Name of stderr error file

module load anaconda
conda activate rare_event
cd /home/x-yliu4/aps/SZ3/python-examples
python benchmark.py -c ./event_detection_data_conf_anvil.yml -m benchmark