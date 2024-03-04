#!/bin/bash
#SBATCH -A ckb.prj
#SBATCH -J pgs_subtype
#SBATCH -o pgs_subtype_%A_%a.out
#SBATCH -e pgs_subtype_%A_%a.err
#SBATCH -p short
#SBATCH -c 1
#SBATCH -a 1-73

line=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' sumstats_files.txt)
sumstats=$(echo $line | cut -f1 -d' ')
build=$(echo $line | cut -f2 -d' ')
neff=$(echo $line | cut -f4 -d' ')
pop=$(echo ${sumstats} | cut -d'.' -f1 | rev | cut -d'_' -f1 | rev | tr '[:upper:]' '[:lower:]')

prscsx=/well/ckb/users/aey472/projects/pgs_subtype/programs/PRScsx/PRScsx.py
prscs=/well/ckb/users/aey472/projects/pgs_subtype/programs/PRScs/PRScs.py

ldscdir=/well/ckb/users/aey472/projects/pgs_subtype/data/sumstats/ldsc

ldref=/well/ckb/users/aey472/projects/pgs_subtype/data/LDref/ldblk_1kg_${pop}

#Rscript fixcolnames.R ${sumstats}
#python harmonise.py ${sumstats} ${build}

bname=$(echo ${sumstats} | cut -d'.' -f1) 

#cd /well/ckb/users/aey472/projects/pgs_subtype/data/sumstats/ldsc
#zgrep -vE 'I|D' ${ldscdir}/${bname}_ldsc.txt.gz | gzip -c > ${ldscdir}/${bname}.ACTG.txt.gz

source /well/ckb/users/aey472/program_files/miniconda3/etc/profile.d/conda.sh
conda activate /gpfs3/well/ckb/users/aey472/program_files/miniconda3/envs/prscsx/

N_THREADS=4
export MKL_NUM_THREADS=${N_THREADS}
export NUMEXPR_NUM_THREADS=${N_THREADS}
export OMP_NUM_THREADS=${N_THREADS}

python ${prscsx} \
  --ref_dir=${ldref} \
  --bim_prefix=VALIDATION_BIM_PREFIX \
  --sst_file=${ldscdir}/${bname}.ACTG.txt.gz \
  --n_gwas=${nef} \
  --pop=${pop} \
  --out_dir=${ldscdir}
  --out_name=${bname}.ACTG
