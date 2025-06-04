This figure was generated using data from Pham, Panda, et al. Cell Stem Cell 2022 https://doi.org/10.1016/j.stem.2022.08.001 .
The data was processed using the scLinaX pipeline from Tomofuji et al. https://github.com/ytomofuji/scLinaX 

We used vartrix instead of cellsnp-lite to pull allelic reads from 10x chromium data, as we found that that works best.

The required input files are position sorted bam files, output from cellranger, and the barcode file, also from cellranger.

You will need annovar, bcftools, and vartrix installed. 
You will also need an R environment with Matrix and Reshape2 libraries

All variables to be changed should be in 1.runVartrix.slurm  . We used a slurm job system, the file will need to be adapted to your system.

After processing with script 1, it will go directly into script 2.
You can then look at our jupyter notebook for downstream processing.