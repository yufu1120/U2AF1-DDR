#!/bin/bash
#SBATCH --job-name=2
#SBATCH -n 6
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=./MV4HDAC8.log
#SBATCH --error=./MV4HDAC8.err

params_dir_ref="$(pwd)"
current_project=./2_MV4HDAC8
current_project_output=./2_MV4HDAC8/output
gtf=./genome.gtf

singularity exec --bind ${params_dir_ref} /opt/rMATS/4.1.2.simg python /rmats/rmats.py \
 --b1 ./lnk/MV4HDAC8_b1_NIR.txt \
 --b2 ./lnk/MV4HDAC8_b2_IR.txt \
 --gtf ${gtf} -t paired --readLength 101 --variable-read-length --nthread 6 \
 --od ${current_project_output}/final_output/ --tmp ${current_project_output}/tmp_output/
