#!/bin/bash
#SBATCH --job-name=homo30_40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH lglassosubject=30_time=40.R asubject=30_time=40.out