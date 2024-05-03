#!/bin/bash

#module load snakemake/5.7.4
#module load r/4.1.1 #load this version on ipop-up, otherwise fails with r/4.0.3

#module load trim-galore/0.6.10
#module load fastqc/0.12.1

#module load star/2.7.5a
##module load bowtie2/2.4.4 #load this version on ipop-up, otherwise fails because lack of libtbb.so.2
#module load samtools/1.13
#module load picard/2.23.5
#module load subread/2.0.1
#module load bedtools/2.30.0

##module load deeptools/3.5.0
##module load htseq/0.13.5
##module load multiqc/1.13


snakemake --unlock
snakemake -nrp -j 200 --resources load=100 --cluster "sbatch -p {cluster.partition} --exclude cpu-node130 --cpus-per-task {cluster.cpu} --mem {cluster.ram} --job-name {rule}" --cluster-config cluster_config.json --latency-wait 1800 --max-jobs-per-second 1 --configfile config.json

