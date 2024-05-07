
## Usage 

Bulk RNAseq data were analysed as follow:

Raw data were placed in a `fastq/raw/` subfolder and processed using the provided snakemake workflows scripts and simple bash script, themselve gathered in a `script/` subfolder. Snakemake workflows scripts are composed of at least 4 files:
- Snakefile: contains the actual workflow script.
- cluster_config.json: contains the resources that are allocated to each workflow's rule. 
- config.json file: must be modified prior to execution to fill in the paths of the working directory, genome indexes and annotation file; adapted to the repository structure of your cluster.
- exeSnakemakeSlurm.sh: bash script to run with the following command line: `sbatch exeSnakemakeSlurm.sh` that launches workflow execution.

These 4 files are accompanied with another required file: 
- genome.txt: a genome file required by bedtools coverage, that defines the expected chromosome order in the input files. 
