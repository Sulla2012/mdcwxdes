#!/bin/bash -l
#SBATCH -t 30:00
#SBATCH -N 8
#SBATCH -q debug
#SBATCH -C haswell
module load anaconda3
source activate mdcwxdes

#cat << EOF > sub_corr_test.sl 
#!/bin/bash
#SBATCH -t 30:00
#SBATCH -N 8
#SBATCH -q debug
#SBATCH -C haswell

python corr_test.py

#EOF

#sbatch sub_corr_test.sl


