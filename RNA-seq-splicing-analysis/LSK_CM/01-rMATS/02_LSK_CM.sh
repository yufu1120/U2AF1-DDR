#!/bin/bash
#SBATCH --job-name=IR
#SBATCH -n 6
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=./LSK_CM.log
#SBATCH --error=./LSK_CM.err

params_dir_ref="$(pwd)"
current_project=.
current_project_output=${current_project}/2
gtf=./refFlat_gene_mm10_200213.gtf


singularity exec --bind ${params_dir_ref} /opt/rMATS/4.1.2.simg python /rmats/rmats.py \
 --b1 ${current_project}/input/lnk/2_CHW_b1_NIR.txt \
 --b2 ${current_project}/input/lnk/2_CHW_b2_IR.txt \
 --gtf ${gtf} -t single --readLength 101 --variable-read-length --nthread 6 --paired-stats \
 --od ${current_project_output}/final_output/ --tmp ${current_project_output}/tmp_output/
