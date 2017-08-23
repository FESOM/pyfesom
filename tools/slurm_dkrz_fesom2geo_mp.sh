#!/bin/bash
#SBATCH --job-name=fesom2geo      # Specify job name
#SBATCH --partition=prepost     # Specify partition name
#SBATCH --ntasks=24            # Specify max. number of tasks to be invoked
#SBATCH --cpus-per-task=2     # Specify number of CPUs per task
#SBATCH --time=10:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=ab0995       # Charge resources on this project account
#SBATCH --output ./my-out.txt
#SBATCH --error ./my-error.txt   # ile name for standard output

# Bind your OpenMP threads
#export OMP_NUM_THREADS=1
#export KMP_AFFINITY=verbose,granularity=core,compact,1
#export KMP_STACKSIZE=64m


python -V
nproc

python fesom2geo_mp.py ./cf/fesom2geo_iexample





