To make this figure, we staretd with the fastq files, our list of known SNPs in H9, and an N-masked genome.

The processing script aligns and deduplicates reads, then splits reads based on which allele they contain and calculates coverage over SNPs.
This script will have to be adapted to your local computing architecture, we use a Slurm based job system so that is what it is written for.

For the plotting script, we did this in R. It should be pretty straightforward.


