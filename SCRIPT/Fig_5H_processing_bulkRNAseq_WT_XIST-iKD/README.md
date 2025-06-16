
## Usage 

Allelic analysis of Bulk RNAseq data was carried out as follow:

Raw data were placed in a `fastq/raw/` subfolder and processed using the provided snakemake workflows scripts and simple bash script, themselve gathered in a `script/` subfolder. The Snakemake workflow is composed of 4 files:
- `Snakefile`: contains the actual workflow script.
- `cluster_config.json`: contains the resources that are allocated to each workflow's rule. 
- `config.json file`: must be modified prior to execution to fill in the paths of the working directory, genome indexes and annotation file; adapted to the repository structure of your cluster.
- `exeSnakemakeSlurm.sh`: bash script to run with the following command line: `sbatch exeSnakemakeSlurm.sh` that launches workflow execution.

This workflow requires other input files: 
- `genome.txt`: a genome file required by bedtools coverage, that defines the expected chromosome order in the input files (given in this folder).
- An N-masked genome to be used for STAR alignment, required by SNPSplit, whose generation instructions is described in 1-Making_SNP-N-masked-genome_for_SNPsplit.md (in this folder), and whose path must be specified in the config.json file prior to execution of the workflow.

The output files `.depthPerSNP` were then processed on R to generate Figure 5H and 5I. 
