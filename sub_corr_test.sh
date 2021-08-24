#!/bin/bash -l

#SBATCH -N 10
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -t 022:30:00

#module load anaconda3
source activate mdcwxdes


srun -n 10 python corr_test.py

#EOF

#sbatch sub_corr_test.sl


