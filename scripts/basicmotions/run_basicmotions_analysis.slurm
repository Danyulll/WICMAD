#!/bin/bash
#SBATCH --job-name=basicmotions-analysis
#SBATCH --account=def-dsteph
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.krasnov@mail.mcgill.ca

# --- Modules: match the versions you used to install packages on the login node ---
module load r/4.5.0

# --- Use the same user library you installed into on the login node ---
export R_LIBS="$HOME/.local/R/$EBVERSIONR/"
export TMPDIR="${SLURM_TMPDIR:-/tmp}"

# --- Keep threaded math libs to 1 thread per PSOCK worker ---
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# (Optional) expose allocation to R for logging
export R_PARALLEL_NUM_WORKERS="${SLURM_CPUS_PER_TASK}"

echo "Job ID: $SLURM_JOB_ID"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "R_LIBS: $R_LIBS"
echo "Starting at: $(date)"
echo "Host: $(hostname)"

# --- Run your simulation ---
srun Rscript run_basicmotions_analysis.R

echo "Finished at: $(date)"
