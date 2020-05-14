FROM ubuntu:latest
MAINTAINER Rob Kitchen <rkitchen@mgh.harvard.edu>

##
## Set paths
##
## Set the path to install dependencies
ARG path_bin=/bin/exceRpt
ENV PATH_BIN=$path_bin




RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc unzip make libbz2-dev zlib1g zlib1g-dev git default-jre libncurses5-dev libncursesw5-dev liblzma-dev libcurl4-openssl-dev


RUN mkdir -p $PATH_BIN \
    && cd $PATH_BIN \
    && echo $PATH_BIN \
    && pwd


## exceRpt
RUN cd $PATH_BIN \
    && git clone https://github.com/rkitchen/exceRpt.git \
    && cp exceRpt/exceRpt_smallRNA . \
    && cp exceRpt/exceRpt_Tools.jar .
ENV EXE_EXCERPT_TOOLS=$path_bin/exceRpt/exceRpt_Tools.jar


## Thunder
RUN cd $PATH_BIN \
    && git clone https://github.com/rkitchen/Thunder.git
ENV EXE_THUNDER=$path_bin/Thunder/Thunder.jar


## STAR
RUN cd $PATH_BIN \
    && git clone https://github.com/alexdobin/STAR.git
ENV EXE_STAR=$path_bin/STAR/bin/Linux_x86_64/STAR

## Salmon
RUN cd $PATH_BIN \
    && wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz \
    && tar xzvf salmon-1.1.0_linux_x86_64.tar.gz
ENV EXE_SALMON=$path_bin/salmon-latest_linux_x86_64/bin/salmon

## Bowtie1
RUN cd $PATH_BIN \
    && wget --no-verbose https://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.2.3/bowtie-1.2.3-linux-x86_64.zip \
    && unzip $PATH_BIN/bowtie-1.2.3-linux-x86_64.zip \
    && rm $PATH_BIN/bowtie-1.2.3-linux-x86_64.zip \
    && mv $PATH_BIN/bowtie-1.2.3-linux-x86_64 $PATH_BIN/bowtie1 
ENV EXE_BOWTIE1=$path_bin/bowtie1/bowtie
ENV EXE_BOWTIE1_BUILD=$path_bin/bowtie1/bowtie-build

## Bowtie2
RUN cd $PATH_BIN \
    && wget --no-verbose https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip \
    && unzip $PATH_BIN/bowtie2-2.3.5.1-linux-x86_64.zip \
    && rm $PATH_BIN/bowtie2-2.3.5.1-linux-x86_64.zip \
    && mv $PATH_BIN/bowtie2-2.3.5.1-linux-x86_64 $PATH_BIN/bowtie2
ENV EXE_BOWTIE2=$path_bin/bowtie2/bowtie2
ENV EXE_BOWTIE2_BUILD=$path_bin/bowtie2/bowtie2-build


## samtools
RUN cd $PATH_BIN \
    && wget --no-verbose https://downloads.sourceforge.net/project/samtools/samtools/1.9/samtools-1.9.tar.bz2 \
    && tar -xvf $PATH_BIN/samtools-1.9.tar.bz2 \
    && rm $PATH_BIN/samtools-1.9.tar.bz2 \
    && cd $PATH_BIN/samtools-1.9 \
    && ./configure \
    && make \
    && cd - \
    && mv $PATH_BIN/samtools-1.9 $PATH_BIN/samtools
ENV EXE_SAMTOOLS=$path_bin/samtools/samtools

##  FastQC
RUN cd $PATH_BIN \
    && wget --no-verbose http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
    && unzip $PATH_BIN/fastqc_v0.11.8.zip \
    && rm $PATH_BIN/fastqc_v0.11.8.zip
ENV EXE_FASTQC=$path_bin/FastQC

## BBduk
RUN cd $PATH_BIN \
    && wget --no-verbose https://downloads.sourceforge.net/project/bbmap/BBMap_38.73.tar.gz \
    && tar -xvf $PATH_BIN/BBMap_38.73.tar.gz \
    && rm $PATH_BIN/BBMap_38.73.tar.gz
ENV PATH_BBMAP=$path_bin/bbmap

## SRAtoolkit
RUN cd $PATH_BIN \
    && wget --no-verbose https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    && tar -xvf $PATH_BIN/sratoolkit.current-ubuntu64.tar.gz \
    && rm $PATH_BIN/sratoolkit.current-ubuntu64.tar.gz
ENV EXE_SRA_TOOLKIT=$path_bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump

## Tidy up
#RUN apt-get remove -y wget gcc bzip2 make zlib1g-dev git \
RUN apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


##
## Add exceRpt executables
##
#ADD exceRpt_smallRNA /exceRpt_bin/exceRpt_smallRNA
#ADD mergePipelineRuns.R /exceRpt_bin/mergePipelineRuns.R
#ADD mergePipelineRuns_functions.R /exceRpt_bin/mergePipelineRuns_functions.R
#ADD exceRpt_Tools.jar /exceRpt_bin/exceRpt_Tools.jar
#ADD LICENSE /exceRpt_bin/LICENSE
#ADD README.md /exceRpt_bin/README.md


##
## Add baseDB and example raw data
##
RUN mv $PATH_BIN/exceRpt/build/DATABASE /exceRpt_DB
RUN mkdir -p /exceRptInput \
    && mv $PATH_BIN/exceRpt/ExampleData/testData_human.fastq.gz /exceRptInput/testData_human.fastq.gz

#ADD exceRpt_coreDB /exceRpt_DB/
#ADD ExampleData/testData_human.fastq.gz /exceRptInput/testData_human.fastq.gz


## To run interactively
#docker run -it --entrypoint /bin/bash rkitchen/excerpt


##
## Entrypoint
##
ENTRYPOINT ["make", "-f", "/bin/exceRpt/exceRpt_smallRNA", "EXE_DIR=/bin/exceRpt", "DATABASE_PATH=/exceRpt_DB", "JAVA_EXE=java", "OUTPUT_DIR=/exceRptOutput", "MAP_EXOGENOUS=off", "N_THREADS=4"]


