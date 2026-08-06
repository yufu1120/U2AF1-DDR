#!/bin/bash
#SBATCH --job-name=HDAC8_OE
#SBATCH -n 6
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=./HDAC8_OE.log
#SBATCH --error=./HDAC8_OE.err

params_dir_ref="$(pwd)"
current_project=./HDAC8_OE
current_project_output=./HDAC8_OE/2
gtf=./refFlat_gene_mm10_200213.gtf

singularity exec --bind ${params_dir_ref} /opt/rMATS/4.1.2.simg python /rmats/rmats.py \
 --b1 ./lnk/HDAC8OE_b1_NIR.txt \
 --b2 ./lnk/HDAC8OE_b2_IR.txt \
 --gtf ${gtf} -t paired --readLength 101 --variable-read-length --nthread 6 \
 --od ${current_project_output}/final_output/ --tmp ${current_project_output}/tmp_output/
