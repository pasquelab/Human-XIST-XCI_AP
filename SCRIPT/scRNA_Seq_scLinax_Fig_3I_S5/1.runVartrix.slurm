#!/bin/bash
#SBATCH --cluster=wice
#SBATCH --job-name="vartrix"
#SBATCH -A lpasquelab
#SBATCH --time=16:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=72
#SBATCH --mem-per-cpu=3400M
#SBATCH --output=%x-%j.log
#SBATCH --get-user-env


#################################




#these variables will need to be edited by you for use:
DIRECTORY=/staging/leuven/stg_00041/Published_StemCell_Data/Pham_EXMC_2022/allelic/
GENOME=/staging/leuven/stg_00041/genomes/GRCh38/GRCh38.p12.genome.fa
BAMfile=timecourse.possorted.bam
BARcodes=barcodes.tsv
VARTRIXloc=/data/leuven/343/vsc34313/vartrix/vartrix_linux
ANNOvarLoc=/staging/leuven/stg_00041/software/annovar/annovar/table_annovar.pl
ANNOVARhumanDB=/staging/leuven/stg_00041/software/annovar/annovar/humandb
RSCRIPTprogram=$VSC_DATA/miniconda3/envs/sclinax/bin/Rscript
RSCRIPTlocation=/staging/leuven/stg_00041/Published_HumanCellAtlas/processvartrixToSclinaxInput.Rscript

#loading modules for our server
module load cluster/wice/batch_icelake
module --force purge
module use /apps/leuven/icelake/2021a/modules/all
module load BCFtools/1.15.1-GCC-10.3.0

cd $DIRECTORY

echo "starting bcftools to find heterozygous SNPs"
bcftools mpileup -f $GENOME -a 'FORMAT/AD' $BAMfile | bcftools filter -O b -i'DP>30' > ${BAM//.bam/.bcf} #calculating variants
bcftools index ${BAM//.bam/.bcf} #indexing bcf file
bcftools call -c -V indels -O b ${BAM//.bam/.bcf} > ${BAM//.bam/.calls.bcf} #calling SNPs
bcftools view -i'QUAL>20' -g het ${BAM//.bam/.calls.bcf} > ${BAM//.bam/.calls.vcf} #filtering for heterozygous SNPs and converting to vcf
bcftools view -H ${BAM//.bam/.calls.vcf} > ${BAM//.bam/.noHeader.vcf}
echo "bcftools complete "

echo "start vartrix "

$VARTRIXloc --bam $BAM \
--cell-barcodes $BARcodes \
--fasta $GENOME \
--out-matrix ${BAM//.bam/.altreads.matrix} \
--out-variants ${BAM//.bam/.variantPos.list} --scoring-method coverage \
--ref-matrix ${BAM//.bam/.refreads.matrix}  \
--vcf ${BAM//.bam/.vartrix.calls.vcf}

echo "vartrix complete"


echo "vartrix complete. starting preprocessing for sclinax input"
date
awk ' { print $1":"$2":"$4":"$5"\t"$2"\t"$4"\t"$5 } '  ${BAM//.bam/.noHeader.vcf} > ${BAM//.bam/.scLinaxfirst4cols.tsv}
awk ' { print $1"\t"$2"\t"$2"\t"$4"\t"$5} ' ${BAM//.bam/.noHeader.vcf} > ${BAM//.bam/.annovarInput.tsv}

echo "Starting Annovar"
$ANNOvarLoc \
${BAM//.bam/.annovarInput.tsv} \
$ANNOVARhumanDB \
-buildver hg38 \
-out ${BAM//.bam/.annovarOutput.tsv} \
-remove -protocol refGene,ensGene,avsnp150,gnomad41_genome \
-operation g,g,f,f -nastring NA -polish

echo "Done Annovar"
date

echo "Starting R script"
#the following R script takes 6 arguments in the order of sclinaxfirst4cols, variantPos, refReads, altReads, barcodes, outfile and sampleName
$RSCRIPTprogram $RSCRIPTlocation ${BAM//.bam/.scLinaxfirst4cols.tsv} ${BAM//.bam/.calls.vcf} ${BAM//.bam/.refreads.matrix} ${BAM//.bam/.testReadGroup.depthPerSNP5.altreads.matrix} $BARcodes ${BAM//.bam/.formattedForscLinax.tsv} ${BAM//.bam/}

echo "script completed"
date