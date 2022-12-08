#!/bin/bash
#SBATCH --job-name=MOLNAME.ATOMNUM.RUNID
#SBATCH --partition=owners,simes
#SBATCH --time=2-00:00:00
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
/home/users/tatang/project/Diamondoid/build/diamondMC
