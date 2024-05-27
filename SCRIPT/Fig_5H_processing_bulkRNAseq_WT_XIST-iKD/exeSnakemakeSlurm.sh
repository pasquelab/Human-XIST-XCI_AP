#!/bin/bash
################################ Slurm options #################################

### Job name
#SBATCH --job-name=Mapping

### Output
#SBATCH --output=Mapping-%j.out  # both STDOUT and STDERR

################################   Modules    #################################

module load snakemake/7.32.4
module list

################################   Workflow    #################################

snakemake --unlock
snakemake -nrp -j 200 --resources load=100 --cluster "sbatch -p {cluster.partition} --exclude cpu-node130 --cpus-per-task {cluster.cpu} --mem {cluster.ram} --job-name {rule}" --cluster-config cluster_config.json --latency-wait 1800 --max-jobs-per-second 1 --configfile config.json --use-singularity

