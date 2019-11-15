#!/bin/sh
### Job name
#PBS -N indexVert
### Number of nodes 
##PBS -l nodes=1:ppn=8
#PBS -l nodes=1:ppn=64
#PBS -q bigmem
#PBS -j oe

cd /gpfs/scratch/fas/gerstein/rrk24/ANNOTATIONS/Genomes_BacteriaFungiMammalPlantProtistVirus

mkdir STAR_GENOME_VERTEBRATE1
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_VERTEBRATE1 --genomeFastaFiles Vertebrate1.fa --limitGenomeGenerateRAM 400000000000 --genomeChrBinNbits 13 --runThreadN 64

mkdir STAR_GENOME_VERTEBRATE2
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_VERTEBRATE2 --genomeFastaFiles Vertebrate2.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 16 --runThreadN 64

mkdir STAR_GENOME_VERTEBRATE3
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_VERTEBRATE3 --genomeFastaFiles Vertebrate3.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 18 --runThreadN 64

mkdir STAR_GENOME_VERTEBRATE4
/gpfs/scratch/fas/gerstein/rrk24/bin/STAR_2.4.0i/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STAR_GENOME_VERTEBRATE4 --genomeFastaFiles Vertebrate4.fa --limitGenomeGenerateRAM 256000000000 --genomeChrBinNbits 18 --runThreadN 64

