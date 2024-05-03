#!/bin/bash
################################ Slurm options #################################

### Job name
#SBATCH --job-name=Downsampling

### Output
#SBATCH --output=Downsampling-%j.out  # both STDOUT and STDERR

################################   Modules    #################################

module load snakemake/7.32.4
module load bedtools/2.30.0 # to be removed when sif updated
module list

################################   Workflow    #################################

snakemake --unlock
snakemake -nrp -j 200 --resources load=100 --cluster "sbatch -p {cluster.partition} --cpus-per-task {cluster.cpu} --mem {cluster.ram} --job-name {rule}" --cluster-config cluster_config.json --latency-wait 1800 --max-jobs-per-second 1 --configfile config.json --use-singularity
# snakemake -n : to perform a dry run 
