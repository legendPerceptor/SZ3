#!/bin/bash

#PBS -N detect_sz3
#PBS -A SDR
#PBS -l walltime=01:00:00
#PBS -l select=1:ngpus=1
#PBS -q gpu
#PBS -o /lcrc/project/ECP-EZ/yuanjian/APS-data/experiment-apr4/detect_sz3.out
#PBS -e /lcrc/project/ECP-EZ/yuanjian/APS-data/experiment-apr4/detect_sz3.err

module load anaconda3/2024.10
conda activate rare_event
cd /home/ac.yuanjian/Research/SZ3/python-examples
python -u detection_benchmark.py