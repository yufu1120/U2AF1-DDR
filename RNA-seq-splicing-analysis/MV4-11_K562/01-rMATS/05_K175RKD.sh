#!/bin/bash
#SBATCH --job-name=5
#SBATCH -n 6
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=./K175RKD.log
#SBATCH --error=./K175RKD.err

params_dir_ref="$(pwd)"
current_project=./5_K175RKD
current_project_output=./5_K175RKD/output
gtf=./genome.gtf

singularity exec --bind ${params_dir_ref} /opt/rMATS/4.1.2.simg python /rmats/rmats.py \
 --b1 ./lnk/K175RKD_b1_NIR.txt \
 --b2 ./lnk/K175RKD_b2_IR.txt \
 --gtf ${gtf} -t paired --readLength 101 --variable-read-length --nthread 6 \
 --od ${current_project_output}/final_output/ --tmp ${current_project_output}/tmp_output/
