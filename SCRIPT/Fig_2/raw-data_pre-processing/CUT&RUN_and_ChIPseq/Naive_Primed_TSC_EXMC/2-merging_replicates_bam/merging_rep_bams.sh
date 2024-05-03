#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=Merging_rep

### Output
#SBATCH --output=Merging_rep-%j.out  # both STDOUT and STDERR

################################   Modules    #################################

#singularity shell cutrun_v0.sif
module load samtools/1.13

input_path=/shared/projects/tsc/cutrun_alldata_test/TSC-EXMC-Primed-Naive/bam/
output_path=/shared/projects/tsc/cutrun_alldata_test/TSC-EXMC-Primed-Naive/bam/merged_downsampled_bam/

if [ ! -d "${output_path}" ]; then
    mkdir -p "${output_path}"
fi


# TSC, EXMC

samtools merge -o ${output_path}TSC_H3K27me3.hg38.dedupped.merged.bam ${input_path}D1539C21.hg38.dedupped.bam ${input_path}D1539C25.hg38.dedupped.bam
samtools merge -o ${output_path}TSC_H3K9me3.hg38.dedupped.merged.bam ${input_path}D1539C22.hg38.dedupped.bam ${input_path}D1539C26.hg38.dedupped.bam
samtools merge -o ${output_path}TSC_H2AK119Ub.hg38.dedupped.merged.bam ${input_path}D1539C23.hg38.dedupped.bam ${input_path}D1539C27.hg38.dedupped.bam
samtools merge -o ${output_path}TSC_IGG.merged.bam ${input_path}D1539C24.hg38.dedupped.bam ${input_path}D1539C28.hg38.dedupped.bam

samtools merge -o ${output_path}EXMC_H3K27me3.hg38.dedupped.merged.bam ${input_path}D1539C29.hg38.dedupped.bam ${input_path}D1539C33.hg38.dedupped.bam
samtools merge -o ${output_path}EXMC_H3K9me3.hg38.dedupped.merged.bam ${input_path}D1539C30.hg38.dedupped.bam ${input_path}D1539C34.hg38.dedupped.bam
samtools merge -o ${output_path}EXMC_H2AK119Ub.hg38.dedupped.merged.bam ${input_path}D1539C31.hg38.dedupped.bam ${input_path}D1539C35.hg38.dedupped.bam
samtools merge -o ${output_path}EXMC_IGG.merged.bam ${input_path}D1539C32.hg38.dedupped.bam ${input_path}D1539C36.hg38.dedupped.bam

# Naive

samtools merge -o ${output_path}PXGL_H3K27me3_Filtered.merged.bam ${input_path}D1249C60_Filtered.bam ${input_path}D1249C76_Filtered.bam
samtools merge -o ${output_path}PXGL_H3K9me3_Filtered.merged.bam ${input_path}D1249C61_Filtered.bam ${input_path}D1249C77_Filtered.bam
samtools merge -o ${output_path}PXGL_H2AK119Ub_Filtered.merged.bam ${input_path}D1249C62_Filtered.bam ${input_path}D1249C78_Filtered.bam
samtools merge -o ${output_path}PXGL_IGG.merged.bam ${input_path}D1249C64_Filtered.bam ${input_path}D1249C80_Filtered.bam

# Primed

samtools merge -o ${output_path}Primed_H3K27me3.hg38.dedupped.merged.bam ${input_path}D725C04.hg38.dedupped.bam ${input_path}D725C11.hg38.dedupped.bam
samtools merge -o ${output_path}Primed_H3K9me3.hg38.dedupped.merged.bam ${input_path}D725C05.hg38.dedupped.bam ${input_path}D725C12.hg38.dedupped.bam
samtools merge -o ${output_path}Primed_H2AK119Ub.hg38.dedupped.merged.bam ${input_path}D725C06.hg38.dedupped.bam ${input_path}D725C13.hg38.dedupped.bam
samtools merge -o ${output_path}Primed_IGG.merged.bam ${input_path}D725C07.hg38.dedupped.bam ${input_path}D725C14.hg38.dedupped.bam



