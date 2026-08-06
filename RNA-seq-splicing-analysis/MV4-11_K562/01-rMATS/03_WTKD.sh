#!/bin/bash
#SBATCH --job-name=3
#SBATCH -n 6
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=./WTKD.log
#SBATCH --error=./WTKD.err

params_dir_ref="$(pwd)"
current_project=./3_WTKD
current_project_output=./3_WTKD/output
gtf=./genome.gtf

singularity exec --bind ${params_dir_ref} /opt/rMATS/4.1.2.simg python /rmats/rmats.py \
 --b1 ./lnk/WTKD_b1_NIR.txt \
 --b2 ./lnk/WTKD_b2_IR.txt \
 --gtf ${gtf} -t paired --readLength 101 --variable-read-length --nthread 6 \
 --od ${current_project_output}/final_output/ --tmp ${current_project_output}/tmp_output/
