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

ldscdir=/well/ckb/users/aey472/projects/pgs_subtype/data/sumstats/ldsc

#Rscript fixcolnames.R ${sumstats}
#python harmonise.py ${sumstats} ${build}

bname=$(echo ${sumstats} | cut -d'.' -f1) 
cd /well/ckb/users/aey472/projects/pgs_subtype/data/sumstats/ldsc
zgrep -vE 'I|D' ${ldscdir}/${bname}_ldsc.txt.gz | gzip -c > ${ldscdir}/${bname}.ACTG.txt.gz
