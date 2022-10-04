#!/bin/bash
#SBATCH --job-name=lglasso20_30_1_80
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH hetersubject=20_time=30_coef=1_nodes=80.R asubject=20_time=30_coef=1_nodes=80.out
