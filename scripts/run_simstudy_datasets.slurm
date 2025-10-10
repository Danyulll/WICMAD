#!/bin/bash
#SBATCH --job-name=simstudy-wicmad
#SBATCH --account=def-dsteph          # e.g., def-abc123
#SBATCH --time=7-00:00:00             # 7 days
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45            # total R workers you want
#SBATCH --mem=160G                    # adjust as needed (e.g., 3–4G per worker)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.krasnov@mail.mcgill.ca

# --- Modules (Alliance clusters) ---
module purge
module load r                         # loads the supported R version

# --- Keep threaded math libs to a single thread per worker (important!) ---
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# --- Scratch/temp + user library (optional but handy) ---
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
export R_LIBS_USER="$HOME/R/${R_VERSION%%.*}/site-library:${R_LIBS_USER}"

echo "Job ID: $SLURM_JOB_ID"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Starting at: $(date)"
echo "Host: $(hostname)"

# Expose allocation to R for logging/override if desired
export R_PARALLEL_NUM_WORKERS="${SLURM_CPUS_PER_TASK}"

# Run your simulation script (the script uses detectCores to set up 45 workers)
srun Rscript run_simstudy_datasets.R

echo "Finished at: $(date)"
