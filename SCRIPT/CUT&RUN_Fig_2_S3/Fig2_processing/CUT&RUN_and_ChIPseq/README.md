## Usage 

CUT&RUN, CUT&Tag and ChIPseq analysis was carried out as follow:

Raw data (fastq files) can be downloaded from GEO with the following accession numbers: 
- Group 1: 
   - CUT&RUN data of TSCs and EXMCs from this paper:
       - [GSM8215698](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215698) TSC_H3K27me3_Cut&Run_Rep1
       - [GSM8215699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215699) TSC_H3K9me3_Cut&Run_Rep1
       - [GSM8215700](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215700) TSC_H2AK199Ub_Cut&Run_Rep1
       - [GSM8215701](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215701) TSC_IGG_Cut&Run_Rep1
       - [GSM8215702](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215702) TSC_H3K27me3_Cut&Run_Rep2
       - [GSM8215703](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215703) TSC_H3K9me3_Cut&Run_Rep2
       - [GSM8215704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215704) TSC_H2AK199Ub_Cut&Run_Rep2
       - [GSM8215705](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215705) TSC_IGG_Cut&Run_Rep2
       - [GSM8215706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215706) EXMC_H3K27me3_Cut&Run_Rep1
       - [GSM8215707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215707) EXMC_H3K9me3_Cut&Run_Rep1
       - [GSM8215708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215708) EXMC_H2AK199Ub_Cut&Run_Rep1
       - [GSM8215709](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215709) EXMC_IGG_Cut&Run_Rep1
       - [GSM8215710](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215710) EXMC_H3K27me3_Cut&Run_Rep2
       - [GSM8215711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215711) EXMC_H3K9me3_Cut&Run_Rep2
       - [GSM8215712](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215712) EXMC_H2AK199Ub_Cut&Run_Rep2
       - [GSM8215713](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8215713) EXMC_IGG_Cut&Run_Rep2
   - CUT&RUN data of Naive and Primed cells from [Alfeghaly et al., 2023](link):
       - [GSM7873828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873828) PXGL, H3K27me3, WT, Rep1
       - [GSM7873829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873829) PXGL, H3K9me3, WT, Rep1
       - [GSM7873830](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873830) PXGL, H2AK119Ub, WT, Rep1
       - [GSM7873831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873831) PXGL, IgG, WT, Rep1
       - [GSM7873840](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873840) PXGL, H3K27me3, WT, Rep2
       - [GSM7873841](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873841) PXGL, H3K9me3, WT, Rep2
       - [GSM7873842](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873842) PXGL, H2AK119Ub, WT, Rep2
       - [GSM7873843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873843) PXGL, IgG, WT, Rep2
       - [GSM7873851](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873851) Primed, H3K27me3, Rep1
       - [GSM7873852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873852) Primed, H3K9me3, Rep1
       - [GSM7873853](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873853) Primed, H2AK119Ub, Rep1
       - [GSM7873854](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873854) Primed, IgG, Rep1
       - [GSM7873855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873855) Primed, H3K4me3, Rep2
       - [GSM7873856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873856) Primed, H3K27me3, Rep2
       - [GSM7873857](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873857) Primed, H3K9me3, Rep2
       - [GSM7873858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873858) Primed, H2AK119Ub, Rep2
       - [GSM7873859](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7873859) Primed, IgG, Rep2
         
- Group 2:
   - ChIPseq data of Naive and Capacitated cells from [Agostinho de Sousa et al., 2023](https://www.science.org/doi/10.1126/sciadv.adg1936): 
       - [GSM6749200](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749200) cR_H9_EOS_d0_rep1, H3K27me3
       - [GSM6749201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749201) cR_H9_EOS_d0_rep2, H3K27me3
       - [GSM6749202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749202) cR_H9_EOS_d10_rep1, H3K27me3
       - [GSM6749203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749203) cR_H9_EOS_d10_rep2, H3K27me3
       - [GSM6749172](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749172) cR_H9_EOS_d0_rep1, H3K9me3
       - [GSM6749173](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749173) cR_H9_EOS_d0_rep2, H3K9me3
       - [GSM6749174](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749174) cR_H9_EOS_d10_rep1, H3K9me3
       - [GSM6749175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6749175) cR_H9_EOS_d10_rep2, H3K9me3
         
- Group 3:
   - CUT&RUN data of H9-NK2 TSCs and CUT&Tag data of CT27 TSCs:
       - [GSM8229337](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8229337) H9 TSC H3K27me3 CUT&RUN
       - [GSM8229338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8229338) CT27 TSC H3K27me3 CUT&Tag

- Group 4:
   - CUT&RUN data of clonal TSCs:
       - [XXXXX](link)
       - [XXXXX](link)


Each group of data was analysed separately in different folder. Fastq files were placed in a `fastq/raw/` subfolder on a cluster and processed using the provided snakemake workflows (as Snakefiles) and simple bash script, in sequential order.
All softwares required for the analysis are shared as a conda environment `.yml` file.
The output files `.regions.bed.gz`, `.idxstats.txt` , `macs2_enlarged.bed` were then processed on R to generate Figure 2B-C-D-E-F-G (see Fig2_R_plotting/ [here](https://github.com/pasquelab/Human-XIST-XCI_AP/tree/main/SCRIPT/Fig_2/Fig2_R_plotting)).