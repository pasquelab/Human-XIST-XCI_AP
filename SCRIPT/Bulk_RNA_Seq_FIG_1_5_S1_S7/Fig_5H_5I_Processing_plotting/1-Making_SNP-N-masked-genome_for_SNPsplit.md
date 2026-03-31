First formatted SNP list into a bed file: 

```
cat H9formattedSNPsplit.txt | awk '{print $2 "\t" $3-1 "\t" $3}' > BBalaton_SNP_H9.bed
```

Then masked SNP in genome fasta file:

```
module load bedtools/2.30.0

maskFastaFromBed -fi /shared/banks/genomes/homo_sapiens/hg38/fasta/hg38.fa -bed BBalaton_SNP_H9.bed -fo N_masked_bbalaton_SNP_hg38.fa
```

Finally generated star index of this N-masked genome:

```
#!/bin/bash

#SBATCH --mem=40G # increase memory to 25G
#SBATCH --cpus-per-task=8

module load star/2.7.5a

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star-index/ --genomeFastaFiles N_masked_bbalaton_SNP_hg38.fa
```
