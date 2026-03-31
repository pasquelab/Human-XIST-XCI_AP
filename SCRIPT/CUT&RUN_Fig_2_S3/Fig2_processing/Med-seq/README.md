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

These fastq files were placed in a `fastq/raw/` subfolder on a cluster and processed using the snakemake workflow provided as a Snakefile in `1_fq_to_bw/`. 
This workflow requires another R script, given in the same folder: `Xeno_human.R`. This file must also be modified prior to execution to set the library paths for R packages (line 3). 

All softwares required for the analysis are shared as a conda environment `.yml` file.

The output med-seq bigwig files were then processed on R to be denoised (see R script given in `2_bw_to_denoised_df_R/`) and the resulting data frames were used to generate Figure 2J (quarto script provided in  `Fig2_R_plotting/` [here](https://github.com/pasquelab/Human-XIST-XCI_AP/tree/main/SCRIPT/Fig_2/Fig2_R_plotting/))

