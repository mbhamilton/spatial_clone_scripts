#!/bin/bash

#SBATCH --job-name="noherb_Circ"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output="noherbx.out"
#SBATCH --error="noherbx.err.out"
#SBATCH --mail-user=cd1231@georgetown.edu
#SBATCH --mail-type=END,FAIL

module load julia

julia noherb.jl