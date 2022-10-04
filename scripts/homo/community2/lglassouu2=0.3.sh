#!/bin/bash
#SBATCH --job-name=homo0.3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH lglassouu2=0.3.R auu2=0.3.out
