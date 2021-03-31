#!/bin/bash
#SBATCH --job-name=MOLNAME_RUNID
#SBATCH --partition=iric,normal
#SBATCH --time=0-02:00:00
##SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
##SBATCH --mem-per-cpu=4000
#SBATCH --error=job.err
#SBATCH --output=job.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=taaatang@gmail.com
# export OMP_NUM_THREADS=4
srun -n 1 /Users/tatang/Documents/work/projects/PPP/build/main.out