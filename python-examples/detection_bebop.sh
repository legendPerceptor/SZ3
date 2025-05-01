#!/bin/bash

#PBS -N detect_sz3
#PBS -A SDR
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=16
#PBS -q bdwall
#PBS -o /lcrc/project/ECP-EZ/yuanjian/APS-data/apr15-metrics/detect_sz3_logscale_A3.out
#PBS -e /lcrc/project/ECP-EZ/yuanjian/APS-data/apr15-metrics/detect_sz3_logscale_A3.err

module load openmpi/4.1.1
module load anaconda3/2024.06
conda activate bebop_sz3
cd /home/ac.yuanjian/Research/SZ3/python-examples
python -u detection_benchmark.py