#!/bin/sh
### Job name
#PBS -N indexPlants
### Number of nodes 
##PBS -l nodes=1:ppn=8
#PBS -l nodes=1:ppn=64
#PBS -q bigmem
#PBS -j oe

cd /gpfs/scratch/fas/gerstein/rrk24/ANNOTATIONS/Genomes_BacteriaFungiMammalPlantProtistVirus

mkdir STAR_GENOME_PLANTS1
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_PLANTS1 --genomeFastaFiles Plant1.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 13 --runThreadN 64

mkdir STAR_GENOME_PLANTS2
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_PLANTS2 --genomeFastaFiles Plant2.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 18 --runThreadN 64

mkdir STAR_GENOME_PLANTS3
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_PLANTS3 --genomeFastaFiles Plant3.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 18 --runThreadN 64

mkdir STAR_GENOME_PLANTS4
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_PLANTS4 --genomeFastaFiles Plant4.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 18 --runThreadN 64

mkdir STAR_GENOME_PLANTS5
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_PLANTS5 --genomeFastaFiles Plant5.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 13 --runThreadN 64
