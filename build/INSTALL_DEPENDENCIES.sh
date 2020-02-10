#!/bin/bash

##
## INSTALL DEPENDENCIES
##

PATH_BIN=~/WORK/exceRpt/build/bin
mkdir -p $PATH_BIN  && cd $PATH_BIN

## exceRpt
git clone https://github.com/rkitchen/exceRpt.git
EXE_EXCERPT_TOOLS=$PATH_BIN/exceRpt/exceRpt_Tools.jar
cp $EXE_EXCERPT_TOOLS $PATH_BIN/
git clone https://github.com/rkitchen/Thunder.git
EXE_THUNDER=$PATH_BIN/Thunder/Thunder.jar

## STAR
git clone https://github.com/alexdobin/STAR.git
EXE_STAR=$PATH_BIN/STAR/bin/MacOSX_x86_64/STAR

## Bowtie1
#wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.2.3/bowtie-1.2.3-linux-x86_64.zip
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.2.3/bowtie-1.2.3-macos-x86_64.zip
unzip $PATH_BIN/bowtie-1.2.3-macos-x86_64.zip
rm $PATH_BIN/bowtie-1.2.3-macos-x86_64.zip
mv $PATH_BIN/bowtie-1.2.3-macos-x86_64 $PATH_BIN/bowtie1
EXE_BOWTIE1=$PATH_BIN/bowtie1/bowtie
EXE_BOWTIE1_BUILD=$PATH_BIN/bowtie1/bowtie-build

## Bowtie2
#wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.5.1/bowtie2-2.3.5.1-macos-x86_64.zip
unzip $PATH_BIN/bowtie2-2.3.5.1-macos-x86_64.zip
rm $PATH_BIN/bowtie2-2.3.5.1-macos-x86_64.zip
mv $PATH_BIN/bowtie2-2.3.5.1-macos-x86_64 $PATH_BIN/bowtie2
EXE_BOWTIE2=$PATH_BIN/bowtie2/bowtie2
EXE_BOWTIE2_BUILD=$PATH_BIN/bowtie2/bowtie2-build

## samtools
wget https://downloads.sourceforge.net/project/samtools/samtools/1.9/samtools-1.9.tar.bz2
tar -xvf $PATH_BIN/samtools-1.9.tar.bz2
rm $PATH_BIN/samtools-1.9.tar.bz2
cd $PATH_BIN/samtools-1.9
./configure
make
cd -
mv $PATH_BIN/samtools-1.9 $PATH_BIN/samtools
EXE_SAMTOOLS=$PATH_BIN/samtools/samtools

##  FastQC
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
unzip $PATH_BIN/fastqc_v0.11.8.zip
rm $PATH_BIN/fastqc_v0.11.8.zip
PATH_FASTQC=$PATH_BIN/FastQC

## BBduk
wget https://downloads.sourceforge.net/project/bbmap/BBMap_38.73.tar.gz
tar -xvf $PATH_BIN/BBMap_38.73.tar.gz
rm $PATH_BIN/BBMap_38.73.tar.gz
PATH_BBMAP=$PATH_BIN/bbmap

## SRAtoolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvf $PATH_BIN/sratoolkit.current-ubuntu64.tar.gz
rm $PATH_BIN/sratoolkit.current-ubuntu64.tar.gz
EXE_SRA_TOOLKIT=$PATH_BIN/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump


