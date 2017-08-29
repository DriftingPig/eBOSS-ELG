#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 20
#SBATCH -t 00:60:00
#SBATCH -L SCRATCH     #note: specify license need for the file systems your job needs, such as SCRATCH,project

srun -n 780 ./NewData_MultiRun.sh 'elg_240_sgc.v2.TSR.SSR.chunk22_subset' 'random-sweep.merged.chunk22_TSR_SSR_subset'
