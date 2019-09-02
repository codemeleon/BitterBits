#!/usr/bin/env bash

# Add relevant channels

conda config --add channels r
conda config --add channels bioconda

# Adding Tool and components

conda install blast blat mafft muscle \
	raxml fasttree \
	bwa bowtie2 mummer kmc\
	virsorter \
	sra-tools spades picard fastqc kraken kraken2 clustalw prokka \
	samtools bcftools unicycler smalt \
	biopython ruffus pysam \
	roary \
	r r-irkernel

conda install -c conda-forge nodejs


pip install recipy ete3 recipy iva
