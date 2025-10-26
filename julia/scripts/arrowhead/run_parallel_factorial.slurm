#!/bin/bash
#SBATCH --job-name=wicmad_factorial
#SBATCH --account=def-dsteph
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-02:00
#SBATCH --output=parallel_factorial_%j.out
#SBATCH --error=parallel_factorial_%j.err

julia parallel_factorial.jl
