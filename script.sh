#!/bin/bash
#SBATCH -A ckb.prj
#SBATCH -J pgs_subtype
#SBATCH -o pgs_subtype_%A_%a.out
#SBATCH -e pgs_subtype_%A_%a.err
#SBATCH -p short
#SBATCH -c 4
#SBATCH -a 1-73

line=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' sumstats_files.txt)
sumstats=$(echo $line | cut -f1 -d' ')
build=$(echo $line | cut -f2 -d' ')

#Rscript fixcolnames.R ${sumstats}

python harmonise.py ${sumstats} ${build}
