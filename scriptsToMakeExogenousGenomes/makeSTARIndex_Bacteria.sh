#!/bin/sh
### Job name
#PBS -N indexBacteria
### Number of nodes 
#PBS -l nodes=1:ppn=64
#PBS -q bigmem
#PBS -j oe

cd /home/rrk24/scratch_grace/ANNOTATIONS/Genomes_BacteriaFungiMammalPlantProtistVirus
STAR="/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0d/STAR"

mkdir STAR_GENOME_BACTERIA1
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA1 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.0 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA2
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA2 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.1 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA3
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA3 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.2 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA4
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA4 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.3 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA5
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA5 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.4 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA6
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA6 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.5 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA7
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA7 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.6 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA8
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA8 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.7 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA9
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA9 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.8 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
mkdir STAR_GENOME_BACTERIA10
$STAR --runMode genomeGenerate --genomeDir STAR_GENOME_BACTERIA10 --genomeFastaFiles BacteriaGenomes_fasta/ensembl_genome.bacteria.fa.9 --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64
