
## Usage 

Allelic analysis of Bulk RNAseq data was carried out as follow:
Raw data were placed in a `fastq/raw/` subfolder and processed using the provided snakemake workflows scripts gathered in a `script/` subfolder.

All softwares required for the analysis are shared as a conda environment `.yml` file.

This workflow requires other input files: 
- `genome.txt`: a genome file required by bedtools coverage, that defines the expected chromosome order in the input files (given in this folder).
- An N-masked genome to be used for STAR alignment, required by SNPSplit, whose generation instructions is described in 1-Making_SNP-N-masked-genome_for_SNPsplit.md (in this folder), and whose path must be specified in the config.json file prior to execution of the workflow.

The output files `.depthPerSNP` were then processed on R to generate Figure 5H and 5I. 
