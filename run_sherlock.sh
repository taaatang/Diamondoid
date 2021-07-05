#!/bin/bash
#SBATCH --job-name=MOLNAMEMOLNUM.RUNID
#SBATCH --partition=iric,normal
#SBATCH --time=0-04:00:00
##SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=4000
#SBATCH --error=job.err
#SBATCH --output=job.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=taaatang@gmail.com
# export OMP_NUM_THREADS=4
srun -n 1 /home/users/tatang/project/Diamondoid/build/main.out
