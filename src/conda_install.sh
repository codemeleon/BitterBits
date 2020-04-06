#!/usr/bin/env bash

# Add relevant channels

conda config --add channels r
conda config --add channels bioconda

# conda install -c anaconda gcc

# Adding Tool and components

conda install blast blat mafft muscle \
	raxml fasttree \
	bwa bowtie2 mummer kmc emboss bedtools \
	virsorter \
	sra-tools spades picard fastqc kraken kraken2 clustalw prokka \
	samtools bcftools unicycler smalt \
	biopython ruffus pysam \
	roary ariba \
	r r-irkernel \
	weblogo

conda install -c conda-forge nodejs


pip install recipy ete3 recipy iva

# conda create -n py365 python=3.6.5
conda activate py365
conda install ariba
conda deactivate
# added /.anmol/anaconda3/envs/py365/bin/ariba in /etc/bash.bashrc


# conda create -n py27 python=2.7.5
conda activate py27
conda install srst2
conda deactivate
# added /.anmol/anaconda3/envs/py27/bin/srst2 in /etc/bash.bashrc
