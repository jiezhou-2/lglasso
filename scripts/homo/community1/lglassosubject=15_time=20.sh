#!/bin/bash
#SBATCH --job-name=homo15_20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH lglassosubject=15_time=20.R asubject=15_time=20.out
