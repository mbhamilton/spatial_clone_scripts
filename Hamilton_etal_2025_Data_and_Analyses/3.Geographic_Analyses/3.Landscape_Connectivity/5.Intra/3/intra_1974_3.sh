#!/bin/bash

#SBATCH --job-name="intra_1974_3_Circ"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output="intra_1974_3x.out"
#SBATCH --error="intra_1974_3x.err.out"
#SBATCH --mail-user=cd1231@georgetown.edu
#SBATCH --mail-type=END,FAIL

module load julia

julia intra_1974_3.jl