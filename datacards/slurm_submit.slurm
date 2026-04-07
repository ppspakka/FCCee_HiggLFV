#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=reixesreis@gmail.com

#SBATCH --qos=cu_htc
#SBATCH --partition=cpugpu
#SBATCH --job-name=Limits
#SBATCH --output=sbatch.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G

python3 -u run_limits.py --threads 32
