#!/bin/bash


################################ Slurm options #################################

### Job name
#SBATCH --job-name=Merging_rep

### Output
#SBATCH --output=Merging_rep-%j.out  # both STDOUT and STDERR


################################   Workflow    #################################

input_path=/shared/projects/tsc/cutrun_alldata_test/Agostinho-de-Sousa-2023_capacitated-data/bam/
output_path=/shared/projects/tsc/cutrun_alldata_test/Agostinho-de-Sousa-2023_capacitated-data/bam/merged_downsampled_bam/


if [ ! -d "${output_path}" ]; then
    mkdir -p "${output_path}"
fi


singularity exec cutrun.sif samtools merge -o ${output_path}H9-d0_H3K27me3_Filtered.merged.bam ${input_path}SRR22367278_Filtered.bam ${input_path}SRR22367277_Filtered.bam
singularity exec cutrun.sif samtools merge -o ${output_path}H9-d0_H3K9me3_Filtered.merged.bam ${input_path}SRR22367282_Filtered.bam ${input_path}SRR22367281_Filtered.bam

singularity exec cutrun.sif samtools merge -o ${output_path}H9-d10_H3K27me3_Filtered.merged.bam ${input_path}SRR22367276_Filtered.bam ${input_path}SRR22367275_Filtered.bam
singularity exec cutrun.sif samtools merge -o ${output_path}H9-d10_H3K9me3_Filtered.merged.bam ${input_path}SRR22367280_Filtered.bam ${input_path}SRR22367279_Filtered.bam


