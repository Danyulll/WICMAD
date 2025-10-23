#!/bin/bash
#SBATCH --job-name=pamap2-analysis
#SBATCH --time=3-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@domain.com

# Load R module
module load R/4.3.0

# Set working directory
cd $SLURM_SUBMIT_DIR

# Run the analysis
srun Rscript run_pamap2_analysis.R
