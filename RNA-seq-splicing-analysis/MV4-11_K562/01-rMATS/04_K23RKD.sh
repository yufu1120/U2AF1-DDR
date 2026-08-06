#!/bin/bash
#SBATCH --job-name=4
#SBATCH -n 6
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=./K23RKD.log
#SBATCH --error=./K23RKD.err

params_dir_ref="$(pwd)"
current_project=./4_K23RKD
current_project_output=./4_K23RKD/output
gtf=./genome.gtf

singularity exec --bind ${params_dir_ref} /opt/rMATS/4.1.2.simg python /rmats/rmats.py \
 --b1 ./lnk/K23RKD_b1_NIR.txt \
 --b2 ./lnk/K23RKD_b2_IR.txt \
 --gtf ${gtf} -t paired --readLength 101 --variable-read-length --nthread 6 \
 --od ${current_project_output}/final_output/ --tmp ${current_project_output}/tmp_output/
