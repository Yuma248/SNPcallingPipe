# SNPcallPipe
Pipeline to call snps using a reference genome

# Download
git clone https://github.com/Yuma248/SNPcallingPipe
  
# Dependecies 

Perl Parallel:::Loops <Enter>

GNU Parallel

Stacks

AdapterRemoval

Bowtie2

BWA

samtools

bcftools

vcftools


The easiest way to install the dependencies is using conda 

conda create --name SNPcallPipe -c conda-forge -c bioconda perl-parallel-loops parallel stacks adapterremoval bowtie2 bwa samtools bcftools vcftools


