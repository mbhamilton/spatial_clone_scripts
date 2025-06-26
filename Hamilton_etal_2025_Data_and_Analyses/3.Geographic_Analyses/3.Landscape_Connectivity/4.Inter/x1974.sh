#!/bin/bash

#SBATCH --job-name="x1974_Circ"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output="x1974x.out"
#SBATCH --error="x1974x.err.out"
#SBATCH --mail-user=cd1231@georgetown.edu
#SBATCH --mail-type=END,FAIL

module load julia

julia x1974.jl