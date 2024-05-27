## Usage 

MeDseq analysis was carried out as follow:

Raw data (fastq files) can be downloaded from GEO with the following accession numbers: 
    - [GSM8215741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215741): Primed MeDseq Rep1
    - [GSM8215742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215742): Primed MeDseq Rep2
    - [GSM8215743](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215743): TSC MeDseq Rep1
    - [GSM8215744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215744): TSC MeDseq Rep2
    - [GSM8215745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215745): EXMC MeDseq Rep1
    - [GSM8215746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215746): EXMC MeDseq Rep2
    - [GSM8215747](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215747): Naive MeDseq Rep1
    - [GSM8215748](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215748): Naive MeDseq Rep2

These fastq files were placed in a `fastq/raw/` subfolder on a cluster and processed using the provided snakemake workflow script. The Snakemake workflow is composed of 4 files:
- `Snakefile`: contains the actual workflow script.
- `cluster_config.json`: contains the resources that are allocated to each workflow's rule. 
- `config.json file`: must be modified prior to execution to fill in the paths of the working directory and genome indexes; adapted to the repository structure of your cluster.
- `exeSnakemakeSlurm.sh`: bash script to run with the following command line: `sbatch exeSnakemakeSlurm.sh` that launches workflow execution.

This workflow requires another input file, given in the same folder:
- Xeno_human.R: an R script called by the workflow. This file must also be modified prior to execution to set the library paths for R packages (line 3). 

All softwares required for the analysis are encapsulated within a Singularity container image provided in XXX, that should be downloaded and placed in the same folder as the snakemake workflow scripts prior to execution.

The output med-seq bigwig files were then processed on R to generate Figure 2J (see link to quarto). All other input files required for plotting are given in `R/input-files/`.  

