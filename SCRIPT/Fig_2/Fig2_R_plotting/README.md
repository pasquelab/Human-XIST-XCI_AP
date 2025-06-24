## Analysis overview

Chromating profiling data (CUT&RUN, Cut&Tag, ChIPseq, MeDseq) analysis was carried out as follow:

Each dataset's raw data (fastq files, available on GEO: [GSE261711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261711) ) were uploaded on a cluster and was firstly pre-processed using snakemake workflows (provided in the `Fig2_processing/` folder [here](https://github.com/pasquelab/Human-XIST-XCI_AP/tree/main/SCRIPT/Fig_2/Fig2_processing/)). 

The resulting output files `.regions.bed.gz`, `.idxstats.txt` , `macs2_enlarged.bed` and med-seq bigwigs were then processed on R to generate Figure 2 panels. All R scripts are given as a quarto document. All other input files required for plotting are given in `R/input-files/`.  

