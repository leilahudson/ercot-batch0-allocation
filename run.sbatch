#!/bin/bash
#SBATCH --job-name=batch0
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G                 # CO on all 126 sites crashed at 64G and at 160G when the MILP started
#SBATCH --time=12:00:00
#SBATCH --output=logs/batch0_%j.out

# submit from the batch0 folder (the job runs where you submit it), after mkdir -p logs
# usage: sbatch run.sbatch                      (full study)
#        MODELS=co3 N_LOADS=50 sbatch run.sbatch  (env vars pass through)

module load julia gurobi
mkdir -p logs results

export SOLVER=gurobi
export THREADS=$SLURM_CPUS_PER_TASK
export RUN_TAG=${RUN_TAG:-texas_${SLURM_JOB_ID}}

julia --project=. run_texas.jl
