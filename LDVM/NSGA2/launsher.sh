#!/bin/bash

                                                                                                                                                                                                                   
#SBATCH -J NSGA2

# SBATCH -N 1
# SBATCH --ntasks-per-node=1

#SBATCH -n 1
#SBATCH --ntasks-per-core=2
#SBATCH --cpus-per-task=4
#SBATCH --time=96:00:00

#SBATCH --partition=long2

#SBATCH --begin=now

#SBATCH --mail-user=brice.martin@isae-supaero.fr
#SBATCH --mail-type=BEGIN,END,FAIL
module load python
source activate pymoo-env
python pareto_optim.py 
