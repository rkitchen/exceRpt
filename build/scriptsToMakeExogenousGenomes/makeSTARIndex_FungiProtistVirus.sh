#!/bin/sh
### Job name
#PBS -N indexFungiProtistVirus
### Number of nodes 
##PBS -l nodes=1:ppn=8
#PBS -l nodes=1:ppn=64
#PBS -q bigmem
#PBS -j oe

cd /home/rrk24/scratch_grace/ANNOTATIONS/Genomes_BacteriaFungiMammalPlantProtistVirus  &&  mkdir STAR_GENOME_FUNGI_PROTIST_VIRUS

/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0d/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_FUNGI_PROTIST_VIRUS --genomeFastaFiles Fungi.fa Protists.fa Virus.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 15 --runThreadN 64
