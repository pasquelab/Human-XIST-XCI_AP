## Analysis overview

Chromating profiling data (CUT&RUN, Cut&Tag, ChIPseq, MeDseq) analysis was carried out as follow:

Each dataset's raw data (fastq, available on GEO: [GSE261711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261711) ) were uploaded on a cluster and was firstly pre-processed using snakemake workflows (provided in the `raw_data_pre-processing/` folder). All softwares required for the analysis are encapsulated within a Singularity container image provided in XXX. The resulting intermediary files were further analyzed using R scripts (available in the `R/` folder). All other input files necessary for the R analysis are provided in `R/input-files/`. 

## Usage 

The different datasets were analysed as groups of data, each group of data was analysed separately in different folder:

    - CUT_RUN_and_ChIPseq/
         - Agostinho-de-Sousa-et-al-2024_data/
         - Naive_Primed_TSC_EXMC/
         - H9-NK2_and_CT27_TSCs/
    - MeDseq/
  
Raw data were placed in a `fastq/raw/` subfolder and processed using the provided snakemake workflows scripts and simple bash script, themselve gathered in a `script/` subfolder. The Snakemake workflows are composed of 4 files:
- `Snakefile`: contains the actual workflow script.
- `cluster_config.json`: contains the resources that are allocated to each workflow's rule. 
- `config.json file`: must be modified prior to execution to fill in the paths of the working directory, genome indexes and annotation file; adapted to the repository structure of your cluster.
- `exeSnakemakeSlurm.sh`: bash script to run with the following command line: `sbatch exeSnakemakeSlurm.sh` that launches workflow execution.

These workflows sometimes require other input files, given in the same folder.  

The output files `.regions.bed.gz`, `.idxstats.txt` , `macs2_enlarged.bed` and med-seq bigwigs were then processed on R to generate Figure 2. All R scripts are given as a quarto document in `R/` folder. All other input files required for plotting are given in `R/input-files/`.  

