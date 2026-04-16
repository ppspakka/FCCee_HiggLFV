#!/bin/bash

#SBATCH --qos=cu_htc
#SBATCH --partition=cpu
#SBATCH --job-name=Limits
#SBATCH --output=sbatch.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G

python3 -u run_limits.py --threads 20
