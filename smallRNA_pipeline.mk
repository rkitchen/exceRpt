#######################################################################################
##                                                                                   ##
##                      ____        _                                                ##
##     _____  _____ ___|  _ \ _ __ | |_                                              ##
##    / _ \ \/ / __/ _ \ |_) | '_ \| __|                                             ##
##   |  __/>  < (_|  __/  _ <| |_) | |_                                              ##
##    \___/_/\_\___\___|_| \_\ .__/ \__|                                             ##
##                           |_|                                                     ##
##                                                                                   ##
##                                                                                   ##
## SmallRNA-seq pipeline - processes a single sequence file from a single sample     ##
##                                                                                   ##
## Author: Rob Kitchen (rob.kitchen@yale.edu)                                        ##
##                                                                                   ##
## Version 3.2.1 (2015-11-02)                                                        ##
##                                                                                   ##
#######################################################################################
EXCERPT_VERSION := 3.2.1


##
## 1) On the command line, be sure to specify the following MANDATORY parameters
##
OUTPUT_DIR              	:= NULL
INPUT_FILE_PATH         	:= NULL
## You can also override the following OPTIONAL parameters on the commandline
CALIBRATOR_LIBRARY      	:= NULL
RANDOM_BARCODE_LENGTH		:= 0
RANDOM_BARCODE_LOCATION 	:= -5p -3p
ENDOGENOUS_LIB_PRIORITY		:= miRNA,tRNA,piRNA,gencode,circRNA


##
## 2) Choose the main organism for smallRNA / genome alignment (hsa + hg19, hsa + hg38, or mmu + mm10)
##
## Human:
#MAIN_ORGANISM_GENOME_ID := hg19
MAIN_ORGANISM_GENOME_ID := hg38
## Mouse:
#MAIN_ORGANISM_GENOME_ID := mm10


##
## 3) Select whether pipeline is run locally, should be 'true' unless this is the Genboree implementation!
##
LOCAL_EXECUTION := true


##
## 4) Choose optional analysis-specific options (or specify these at the command line)
##
ADAPTER_SEQ                     := guessKnown
REMOVE_LARGE_INTERMEDIATE_FILES := false
SAMPLE_NAME             		:= NULL
QFILTER_MIN_READ_FRAC           := 80
QFILTER_MIN_QUAL                := 20



##
## 5) Choose what kind of EXOGENOUS alignments to attempt 
## 		- off 		: none
##		- miRNAs	: map only to exogenous miRNAs in miRbase
##		- on		: map to exogenous miRNAs in miRbase AND the genomes of all sequenced species in ensembl/NCBI
#MAP_EXOGENOUS      := off
#MAP_EXOGENOUS      := miRNA
MAP_EXOGENOUS      := on


##
### If this is a local installation of the pipeline, be sure to also modify the parameters in steps 4, 5, and 6 below...
##
ifeq ($(LOCAL_EXECUTION),true)

	##
	## 5) Modify installation-specific variables
	##
	N_THREADS 			:= 8
	JAVA_RAM  			:= 64G
	MAX_RAM   			:= 64000000000
	BOWTIE_CHUNKMBS 	:= 8000
	SAMTOOLS_SORT_MEM 	:= 2G
	## NB: The 'EXE_DIR' MUST be an ABSOLUTE PATH or sRNABench will fail!
	EXE_DIR   := /gpfs/scratch/fas/gerstein/rrk24/bin/smallRNAPipeline
	
	
	##
	## 6) Check that the paths to the required 3rd party executables work!
	##
	JAVA_EXE              := /usr/bin/java
	FASTX_CLIP_EXE        := $(EXE_DIR)/fastx_0.0.14/bin/fastx_clipper
	FASTX_FILTER_EXE      := $(EXE_DIR)/fastx_0.0.14/bin/fastq_quality_filter
	BOWTIE1_EXE           := $(EXE_DIR)/bowtie-1.1.1/bowtie
	#VIENNA_PATH          := $(EXE_DIR)/ViennaRNA_2.1.5/bin
	BOWTIE2_EXE           := $(EXE_DIR)/bowtie2-2.2.4/bowtie2
	SAMTOOLS_EXE          := $(EXE_DIR)/samtools-1.1/samtools
	FASTQC_EXE            := $(JAVA_EXE) -classpath $(EXE_DIR)/FastQC_0.11.2:$(EXE_DIR)/FastQC_0.11.2/sam-1.103.jar:$(EXE_DIR)/FastQC_0.11.2/jbzip2-0.9.jar
	SRATOOLS_EXE          := $(EXE_DIR)/sratoolkit.2.5.1-centos_linux64/bin/fastq-dump
	THUNDER_EXE           := $(EXE_DIR)/Thunder.jar
	DATABASE_PATH         := $(EXE_DIR)/DATABASE
	STAR_EXE              := $(EXE_DIR)/STAR_2.4.2a/bin/Linux_x86_64/STAR
	STAR_GENOMES_DIR      := /gpfs/scratch/fas/gerstein/rrk24/ANNOTATIONS/Genomes_BacteriaFungiMammalPlantProtistVirus
	STAR_PARAMS_FILE_PATH := $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in

	##
	## Use the input path to infer filetype and short name
	##
	INPUT_FILE_NAME := $(notdir $(INPUT_FILE_PATH))
	INPUT_FILE_ID   := $(basename $(INPUT_FILE_NAME))
	
else
	##
	## These parameters are for the Genboree installation only
	##
	EXE_DIR 			:= $(SCRATCH_DIR)
	N_THREADS 			:= $(N_THREADS)
	JAVA_RAM 			:= 64G
	MAX_RAM  			:= 64000000000
	BOWTIE_CHUNKMBS 	:= 8000
	SAMTOOLS_SORT_MEM 	:= 2G

	FASTX_CLIP_EXE := fastx_clipper
	FASTX_FILTER_EXE := fastq_quality_filter
	BOWTIE1_EXE := bowtie
	VIENNA_PATH := NULL
	BOWTIE2_EXE := bowtie2
	SAMTOOLS_EXE := samtools
	FASTQC_EXE := $(JAVA_EXE) -classpath $(FASTQC_EXE_DIR):$(FASTQC_EXE_DIR)/sam-1.103.jar:$(FASTQC_EXE_DIR)/jbzip2-0.9.jar
	SRATOOLS_EXE := fastq-dump
	SRNABENCH_EXE := $(SRNABENCH_EXE)
	THUNDER_EXE := $(THUNDER_EXE)
	DATABASE_PATH := $(EXCERPT_DATABASE)
	
	## Path to sRNABench libraries
	SRNABENCH_LIBS := $(SRNABENCH_LIBS)

	STAR_EXE := STAR
	STAR_GENOMES_DIR := $(STAR_GENOMES_DIR)
	STAR_PARAMS_FILE_PATH := $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	
	INPUT_FILE_NAME := $(notdir $(INPUT_FILE_PATH))
    INPUT_FILE_ID := $(INPUT_FILE_ID)
endif



## Define current time
ts := `/bin/date "+%Y-%m-%d--%H:%M:%S"`
## Define tool name for reporting
SMRNAPIPELINE := exceRpt_smallRNA


##
## Initialise smallRNA alignment parameters
##
BOWTIE_SEED_LENGTH := 19

## for endogenous alignments (can be 0-3)
BOWTIE1_MAX_MISMATCHES := 1

## for calibrator, contaminant, rRNA, and exogenous alignments (can only be 0 or 1)
BOWTIE2_MAX_MISMATCHES := 1

# Max number of multiaps to report when aligning to the genome (large value could result in the alignment taking a very long time!)
N_GENOME_MULTIMAPS_TO_REPORT := 10

STAR_outFilterMismatchNoverLmax := 0.05
#BOWTIE_LIB_PATH := $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)





USEAGE := 
ifeq ($(INPUT_FILE_ID),NULL)
  #USEAGE := "make -f smallRNA_pipeline INPUT_FILE_PATH=[required: absolute/path/to/input/.fa|.fq|.sra] N_THREADS=[required: number of threads] OUTPUT_DIR=<required: absolute/path/to/output> INPUT_FILE_ID=[required: samplename] ADAPTER_SEQ=[optional: will guess sequence if not provided here; none, if already clipped input] MAIN_ORGANISM=[optional: defaults to 'hsa'] MAIN_ORGANISM_GENOME_ID=[optional: defaults to 'hg38'] CALIBRATOR_LIBRARY=[optional: path/to/bowtie/index/containing/calibrator/sequences] TRNA_MAPPING=[optional: TRUE|FALSE, default is TRUE] GENCODE_MAPPING=[optional: TRUE|FALSE, default is TRUE] PIRNA_MAPPING=[optional: TRUE|FALSE, default is TRUE] MAP_EXOGENOUS=[optional: off|miRNA|on, default is miRNA]"
  USEAGE := "make -f smallRNA_pipeline INPUT_FILE_PATH=[required: absolute/path/to/input/.fa|.fq|.sra] N_THREADS=[required: number of threads] OUTPUT_DIR=<required: absolute/path/to/output> INPUT_FILE_ID=[required: samplename] ADAPTER_SEQ=[optional: will guess sequence if not provided here; none, if already clipped input] MAIN_ORGANISM=[optional: defaults to 'hsa'] MAIN_ORGANISM_GENOME_ID=[optional: defaults to 'hg38'] CALIBRATOR_LIBRARY=[optional: path/to/bowtie/index/containing/calibrator/sequences] MAP_EXOGENOUS=[optional: off|miRNA|on, default is miRNA]"
endif




##
## Map reads to plant and virus miRNAs
##
ifeq ($(MAP_EXOGENOUS),miRNA)		## ALIGNMENT TO ONLY EXOGENOUS MIRNA
	PROCESS_SAMPLE_REQFILE := EXOGENOUS_miRNA/unaligned.fq.gz
	#PROCESS_SAMPLE_REQFILE := EXOGENOUS_rRNA/unaligned.fq.gz
else ifeq ($(MAP_EXOGENOUS),on)	## COMPLETE EXOGENOUS GENOME ALIGNMENT
	PROCESS_SAMPLE_REQFILE := EXOGENOUS_genomes/ExogenousGenomicAlignments.result.txt
else
	PROCESS_SAMPLE_REQFILE := endogenousUnaligned_ungapped_noLibs.fq
endif


##
## List of plant and virus species IDs to which to map reads that do not map to the genome of the primary organism
##
#EXOGENOUS_MIRNA_SPECIES := $(shell cat $(SRNABENCH_LIBS)/libs/mature.fa | grep ">" | awk -F '-' '{print $$1}' | sed 's/>//g'| sort | uniq | tr '\n' ':' | rev | cut -c 2- | rev)

## Parameters to use for the bowtie mapping of calibrator oligos and rRNAs
BOWTIE2_MAPPING_PARAMS_CALIBRATOR := -D 15 -R 2 -N $(BOWTIE2_MAX_MISMATCHES) -L $(BOWTIE_SEED_LENGTH) -i S,1,0
BOWTIE2_MAPPING_PARAMS_RRNA       := -D 15 -R 2 -N $(BOWTIE2_MAX_MISMATCHES) -L $(BOWTIE_SEED_LENGTH) -i S,1,0


#################################################





##
## Generate unique ID from the input fastq filename and user's sample ID
##
SAMPLE_ID := $(INPUT_FILE_ID)
ifneq ($(SAMPLE_NAME),NULL)
  SAMPLE_ID := $(SAMPLE_ID)_$(SAMPLE_NAME)
endif



##
## Detect filetype and extract from SRA format if necessary
##
COMMAND_CONVERT_SRA := cat $(INPUT_FILE_PATH)
ifeq ($(suffix $(INPUT_FILE_NAME)),.sra)
	COMMAND_CONVERT_SRA := $(SRATOOLS_EXE) --stdout $(INPUT_FILE_PATH)
else ifeq ($(suffix $(INPUT_FILE_NAME)),.gz)
	COMMAND_CONVERT_SRA := gunzip -c $(INPUT_FILE_PATH)
endif


##
## Guess quality encoding
##
#COMMAND_FILTER_BY_QUALITY ?= gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz | $(FASTX_FILTER_EXE) -v -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.tmp 2>>$(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Logic block to write the adapter sequence (whether or not one is provided by the user) to the .adapterSeq file
##
ifeq ($(ADAPTER_SEQ),NULL)
	COMMAND_WRITE_ADAPTER_SEQ := $(COMMAND_CONVERT_SRA) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 40000000 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | $(BOWTIE2_EXE) --no-head -p $(N_THREADS) --local -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 -k 2 --upto 10000000 -x $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/$(MAIN_ORGANISM_GENOME_ID) -U - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if ($$5==255) print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.unique.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.unique.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{print $$6}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort -rnk 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 100 > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).cigarFreqs 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.unique.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{if ($$2==0) print $$3"\t"$$4"\t"$$6"\t"$$10}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | grep "[[:space:]]2[0-9]M[0-9][0-9]S" > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{print $$3}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort -rnk 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 100 > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).okCigarFreqs 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).okCigarFreqs 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{print substr($$2,1,2)}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.txt 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | grep "[[:space:]]$$(<$(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.txt)" 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{getline len<"$(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.txt"; print substr($$4,len+1)}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sed 's/[A]*$$//' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort -rnk 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{if ($$1 > 75) print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).potentialAdapters.txt 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	head -n 1 $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).potentialAdapters.txt | awk '{print $$2}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq;  \
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.*
	LOGENTRY_WRITE_ADAPTER := $(ts) $(SMRNAPIPELINE): Identifying unknown 3' adapter sequence. Removing 3' adapter sequence using fastX:\n
else ifeq ($(ADAPTER_SEQ),guessKnown)
	COMMAND_WRITE_ADAPTER_SEQ := $(COMMAND_CONVERT_SRA) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) FindAdapter -n 10000 -m 1000000 -s 7 -a $(DATABASE_PATH)/adapters/adapters.fa - > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	LOGENTRY_WRITE_ADAPTER := $(ts) $(SMRNAPIPELINE): Identifying 3' adapter from list of known sequences.  Removing 3' adapter sequence using fastX:\n
else ifeq ($(ADAPTER_SEQ),none)
	COMMAND_WRITE_ADAPTER_SEQ := echo 'no adapter' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq;
	COMMAND_CLIP_ADAPTER := $(COMMAND_CONVERT_SRA) | gzip -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.tmp.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.tmp.gz | wc -l | awk '{print "input\t"$$0/4"\nsuccessfully_clipped\tNA"}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	LOGENTRY_WRITE_ADAPTER := Provided 3' adapter clipped input sequence file. No clipping necessary.\n 
else
	COMMAND_WRITE_ADAPTER_SEQ := echo $(ADAPTER_SEQ) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq
	LOGENTRY_WRITE_ADAPTER := $(ts) $(SMRNAPIPELINE): Provided 3' adapter sequence. Removing 3' adapter sequence using fastX:\n
endif


## If no adapter clipping command has been set- use this one:
COMMAND_CLIP_ADAPTER ?= $(COMMAND_CONVERT_SRA) | $(FASTX_CLIP_EXE) -a $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq) -l 15 -v -n -M 7 -z -o $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.tmp.gz >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err



ifeq ($(RANDOM_BARCODE_LENGTH),0)
	COMMAND_REMOVE_RANDOM_BARCODE := mv $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.tmp.gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz
	ENDOGENOUS_QUANT_RANDOM_BARCODE_STATS := 
else
	#COMMAND_REMOVE_RANDOM_BARCODE := gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.tmp.gz | $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) ProcessFastqWithRandomBarcode -n $(RANDOM_BARCODE_LENGTH) $(RANDOM_BARCODE_LOCATION) -stats $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.barcodeStats - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | gzip -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz
	#ENDOGENOUS_QUANT_RANDOM_BARCODE_STATS := -randombarcode $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.barcodeStats
	COMMAND_REMOVE_RANDOM_BARCODE := gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.tmp.gz | $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) ProcessFastqWithRandomBarcode -n $(RANDOM_BARCODE_LENGTH) $(RANDOM_BARCODE_LOCATION) - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | gzip -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz
	ENDOGENOUS_QUANT_RANDOM_BARCODE_STATS := 
endif


##
## Logic block for removing rRNAs and [optionally] calibrator sequences that may have been spiked into the sample
##
ifeq ($(CALIBRATOR_LIBRARY),NULL)
	
	LOGENTRY_MAP_CALIBRATOR_1 := No calibrator sequences\n 
	LOGENTRY_MAP_CALIBRATOR_2 := Moving on to UniVec and rRNA sequences\n
	COMMAND_COUNT_CALIBRATOR := echo -e "calibrator\tNA" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	COMMAND_MAP_CALIBRATOR := 
	
	FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT := $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz
	
else
	
	COMMAND_COUNT_CALIBRATOR := cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.calibratormapped.counts | awk '{sum+=$$1} END {print "calibrator\t"sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	COMMAND_MAP_CALIBRATOR := $(BOWTIE2_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_CALIBRATOR) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noCalibrator.fastq.gz -x $(CALIBRATOR_LIBRARY) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | tee $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.calibratormapped.bam | $(SAMTOOLS_EXE) view - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.calibratormapped.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	LOGENTRY_MAP_CALIBRATOR_1 := $(ts) $(SMRNAPIPELINE): Mapping reads to calibrator sequences using bowtie:\n
	LOGENTRY_MAP_CALIBRATOR_2 := $(ts) $(SMRNAPIPELINE): Finished mapping to the calibrators\n
	
	FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT := $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noCalibrator.fastq.gz
	
endif


##
## Bowtie2 command to align reads to the UniVec contaminant sequence database
##
COMMAND_MAP_UNIVEC := $(BOWTIE2_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noUniVecContaminants.fastq.gz -x $(DATABASE_PATH)/UniVec/UniVec_Core.contaminants -U $(FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam | grep -v "^@" | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminants.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam | grep -v "^@" | awk '{print $$1}' | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c | wc -l > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminants.readCount 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.bam; \
rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam


##
## Bowtie2 command to align reads to the rRNA sequences
##
COMMAND_MAP_RRNAS := $(BOWTIE2_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz -x $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/$(MAIN_ORGANISM_GENOME_ID)_rRNA -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noUniVecContaminants.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam | grep -v "^@" | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNA.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam | grep -v "^@" | awk '{print $$1}' | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c | wc -l > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNA.readCount 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.bam; \
rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam



##
## Remove some potentially large intermediate pipeline output (can save as much as 50% total output size)
##
TIDYUP_COMMAND := 
ifeq ($(REMOVE_LARGE_INTERMEDIATE_FILES),true)
	TIDYUP_COMMAND := rm $(OUTPUT_DIR)/$(SAMPLE_ID)/genome.parsed; rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped*.fastq.gz
endif


##
## Compress only the most vital output!
##
#ls -lh $(OUTPUT_DIR)/$(SAMPLE_ID)/noGenome | awk '{print $$9}' | grep "sense.grouped\|stat" | awk '{print "$(SAMPLE_ID)/noGenome/"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt;
COMPRESS_COMMAND := ls -lh $(OUTPUT_DIR)/$(SAMPLE_ID) | awk '{print $$9}' | grep "readCounts_\|.readLengths.txt\|_fastqc.zip\|.counts\|.adapterSeq\|.qualityEncoding" | awk '{print "$(SAMPLE_ID)/"$$1}' > $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt; \
echo $(SAMPLE_ID).log >> $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt; \
echo $(SAMPLE_ID).stats >> $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt; \
ls -lh $(OUTPUT_DIR)/$(SAMPLE_ID) | awk '{print $$9}' | grep "calibratormapped.counts" | awk '{print "$(SAMPLE_ID)/"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt; \
ls -lh $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes | awk '{print $$9}' | grep "ExogenousGenomicAlignments.result.txt" | awk '{print "$(SAMPLE_ID)/EXOGENOUS_genomes/"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt

#ls -lh $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA | awk '{print $$9}' | grep "exogenous_miRBase_mapped" | awk '{print "$(SAMPLE_ID)/EXOGENOUS_miRNA/"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt
#ls -lh $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA | awk '{print $$9}' | grep "exogenous_miRBase_mapped" | awk '{print "$(SAMPLE_ID)/EXOGENOUS_rRNA/"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt




###########################################################
###########################################################
###########################################################

##
## Main make target
##
.PHONY: all
.DEFAULT: all
all: processSample


##
## Target to selectively compress only the most useful results for downstream processing
##
## - this will typically reduce the volume of data needing to be transferred by 100x
##
compressCoreResults:
	$(COMPRESS_COMMAND)
	tar -cvz -C $(OUTPUT_DIR) -T $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt -f $(OUTPUT_DIR)/$(SAMPLE_ID)_CORE_RESULTS_v$(EXCERPT_VERSION).tgz 2> /dev/null
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt


##
## Delete sample results and logfiles
##
#clean: 
#	rm -r $(OUTPUT_DIR)/$(SAMPLE_ID)




####
#### Main sub-target
####
processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(PROCESS_SAMPLE_REQFILE)
	## Wrap up logging and stats files
	@echo -e "$(ts) $(SMRNAPIPELINE): END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): END\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "#END OF STATS from the exceRpt smallRNA-seq pipeline. Run completed at $(ts)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Compress core results files automatically
	$(COMPRESS_COMMAND)
	#tar -cvz -C $(OUTPUT_DIR) -T $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt -f $(OUTPUT_DIR)/$(SAMPLE_ID)_results.tgz 2> /dev/null
	tar -cvz -C $(OUTPUT_DIR) -T $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt -f $(OUTPUT_DIR)/$(SAMPLE_ID)_CORE_RESULTS_v$(EXCERPT_VERSION).tgz 2> /dev/null
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt
	## END PIPELINE




##
###
####  BEGIN PIPELINE
####
####  vvv Sub-targets to do the read-preprocessing, calibrator mapping, rRNA mapping, en-exRNA mapping, and ex-exRNA mapping vvv
###
##


##
## Make results directory & Write adapter sequence
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat: 
	#$(EXPORT_CMD)
	@echo -e "$(USEAGE)"
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "$(ts) $(SMRNAPIPELINE): BEGIN exceRpt smallRNA-seq pipeline v.$(EXCERPT_VERSION) for sample $(SAMPLE_ID)\n======================\n" > $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): BEGIN \n" > $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) $(SMRNAPIPELINE): Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Processing adapter sequence:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_WRITE_ADAPTER_SEQ)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_WRITE_ADAPTER_SEQ)
	@echo -e "$(ts) $(SMRNAPIPELINE): Progress_1_FoundAdapter" > $(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat
	#
	@echo -e "#STATS from the exceRpt smallRNA-seq pipeline v.$(EXCERPT_VERSION) for sample $(SAMPLE_ID). Run started at $(ts)" > $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	@echo -e "Stage\tReadCount" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## CLIP 3' adapter sequence
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat
	## Run the SW alignment of known adapters regardless of user preference
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Checking adapter against known sequences: $(COMMAND_CONVERT_SRA) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) FindAdapter -n 1000 -m 100000 -s 7 -a $(DATABASE_PATH)/adapters/adapters.fa - > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).knownAdapterSeq 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_CONVERT_SRA) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) FindAdapter -n 1000 -m 100000 -s 7 -a $(DATABASE_PATH)/adapters/adapters.fa - > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).knownAdapterSeq 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Known adapter sequence: $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).knownAdapterSeq)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Carry on with the adapter provided / guessed
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Provided/guessed adapter sequence: $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(LOGENTRY_WRITE_ADAPTER)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_CLIP_ADAPTER)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_CLIP_ADAPTER)
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished removing adapters\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count reads input to adapter clipping
	grep "Input: " $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "input\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Count reads output following adapter clipping
	grep "Output: " $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "successfully_clipped\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Remove random barcodes if there are any
	@echo -e "$(ts) $(SMRNAPIPELINE): Removing $(RANDOM_BARCODE_LENGTH)N random barcode:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_REMOVE_RANDOM_BARCODE)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_REMOVE_RANDOM_BARCODE)
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished removing random barcode\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Guess Fastq quality encoding
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Guessing encoding of fastq read-qualities:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_CONVERT_SRA) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#$(COMMAND_CONVERT_SRA) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@
	@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_CONVERT_SRA) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | grep -v "^\*" | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=75) print "33"; else if(max>74 && min>=64) print "64"; else if(min>=59 && min<64 && max>74) print "64"; else print "64";}' > $@\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_CONVERT_SRA) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | grep -v "^\*" | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=75) print "33"; else if(max>75 && min>=64) print "64"; else if(min>=59 && min<64 && max>74) print "64"; else print "64";}' > $@
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished guessing encoding of fastq read-qualities:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## FILTER clipped reads that have poor overall base quality  &  Remove homopolymer repeats
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Filtering reads by base quality:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_FILTER_BY_QUALITY)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	#$(COMMAND_FILTER_BY_QUALITY)
	@echo -e "$(ts) $(SMRNAPIPELINE): gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz | $(FASTX_FILTER_EXE) -v -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.tmp\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz | $(FASTX_FILTER_EXE) -v -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.tmp 2>>$(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished filtering reads by base quality\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count reads that failed the quality filter
	grep "low-quality reads" $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "failed_quality_filter\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	#
	# Filter homopolymer reads (those that have too many single nt repeats)
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Filtering homopolymer repeat reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) RemoveHomopolymerRepeats -m 0.66 -i $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.tmp -o $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) RemoveHomopolymerRepeats --verbose -m 0.66 -i $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.tmp -o $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.REMOVEDRepeatReads.fastq
	gzip $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished filtering homopolymer repeat reads\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count homopolymer repeat reads that failed the quality filter
	grep "Done.  Sequences removed" $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk -F "=" '{print "failed_homopolymer_filter\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## Assess Read-lengths after clipping
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.readLengths.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz
	@echo -e "======================" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Calculating length distribution of clipped reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) GetSequenceLengths $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq > $@ 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	gunzip -c $< > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq
	$(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) GetSequenceLengths $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq > $@ 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished calculating read-lengths\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Perform FastQC after adapter removal
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered_fastqc.zip: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz
	@echo -e "======================" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Running FastQC on clipped reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/bbsrc/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished running FastQC\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## MAP to external bowtie (calibrator?) library and to UniVec sequences
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noUniVecContaminants.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.readLengths.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered_fastqc.zip
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(LOGENTRY_MAP_CALIBRATOR_1)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_MAP_CALIBRATOR)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_CALIBRATOR)
	@echo -e "$(ts) $(SMRNAPIPELINE): $(LOGENTRY_MAP_CALIBRATOR_2)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count calibrator oligo reads
	$(COMMAND_COUNT_CALIBRATOR)
	#
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to contaminant sequences in UniVec using Bowtie2:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_MAP_UNIVEC)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_UNIVEC)
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to the UniVec contaminant DB\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count UniVec contaminant reads
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminants.readCount | awk '{print "UniVec_contaminants\t"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	

##
## MAP to rRNA sequences
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noUniVecContaminants.fastq.gz
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to ribosomal RNA sequences using Bowtie2:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(COMMAND_MAP_RRNAS)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_RRNAS) 
	$(SAMTOOLS_EXE) sort $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sorted
	$(SAMTOOLS_EXE) index $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sorted.bam
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.bam
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to the rRNAs\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count rRNA reads
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNA.readCount | awk ' {print "rRNA\t"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## Perform FastQC again after rRNA / UniVec removal
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA_fastqc.zip: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz
	@echo -e "======================" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Running FastQC on cleaned reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/bbsrc/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished running FastQC\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log



##
## Map reads to the endogenous genome and transcriptome
##

## map ALL READS to the GENOME (bowtie 1 ungapped)
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA_fastqc.zip
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to the genome of the primary organism ($(MAIN_ORGANISM_GENOME_ID)):\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) -k $(N_GENOME_MULTIMAPS_TO_REPORT) --sam --fullref --best --strata -e 2000 --al $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq --un $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/$(MAIN_ORGANISM_GENOME_ID) - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.bowtie1stats | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) -k $(N_GENOME_MULTIMAPS_TO_REPORT) --sam --fullref --best --strata -e 2000 --al $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq --un $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/$(MAIN_ORGANISM_GENOME_ID) - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.bowtie1stats | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count reads input to genome and transcriptome alignment
	grep "# reads processed:" $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.bowtie1stats | awk -F ": " '{print "reads_used_for_alignment\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Count reads mapped to the genome
	grep "# reads with at least one reported alignment:" $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.bowtie1stats | awk -F ": " '{print $$2}' | awk -F " " '{print "genome\t"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to the genome of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## sort genomic alignments
##$(OUTPUT_DIR)/$(SAMPLE_ID)/balls: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
#$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.unsorted.bam
#	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
#	@echo -e "$(ts) $(SMRNAPIPELINE): Sorting genomic alignments:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
#	@echo -e "$(ts) $(SMRNAPIPELINE): $(SAMTOOLS_EXE) sort -m 8G -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.unsorted.bam > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
#	$(SAMTOOLS_EXE) sort -m 8G -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.unsorted.bam > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
#	@echo -e "$(ts) $(SMRNAPIPELINE): Finished sorting genomic alignments\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## convert genomic alignments to wiggle file for display - CURRENTLY DISABLED - REQUIRES BAM TO BE SORTED!
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.wig: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Converting genomic alignments to .wig format:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(SAMTOOLS_EXE) view -h $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam | $(SAMTOOLS_EXE) mpileup -t SP - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | perl -ne 'BEGIN{print "track type=wiggle_0 name=$(SAMPLE_ID) description=$(SAMPLE_ID)\n"};($$c, $$start, undef, $$depth) = split; if ($c ne $$lastC) { print "variableStep chrom=$$c\n"; };$$lastC=$$c;next unless $$. % 10 ==0;print "$$start\t$$depth\n" unless $$depth<3;' > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.wig 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	$(SAMTOOLS_EXE) view -h $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam | $(SAMTOOLS_EXE) mpileup -t SP - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | perl -ne 'BEGIN{print "track type=wiggle_0 name=$(SAMPLE_ID) description=$(SAMPLE_ID)\n"};($$c, $$start, undef, $$depth) = split; if ($$c ne $$lastC) { print "variableStep chrom=$$c\n"; };$$lastC=$$c;next unless $$. % 10 ==0;print "$$start\t$$depth\n" unless $$depth<3;' > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.wig 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished converting genomic alignments to .wig format\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## map ALL READS to the PRECURSORs (bowtie 1 ungapped)
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_miRNA_precursor.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to miRNA precursors of the primary organism:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that map to the genome to the miRNA library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRNA_precursors - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_miRNA_precursor.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRNA_precursors - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_miRNA_precursor.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that DO NOT map to the genome to the miRNA library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRNA_precursors - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_miRNA_precursor.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRNA_precursors - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_miRNA_precursor.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to miRNA precursors of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## map ALL READS to tRNAs
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_tRNA.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to tRNAs of the primary organism:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that map to the genome to the tRNA library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/tRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_tRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/tRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - 2> /dev/null | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_tRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that DO NOT map to the genome to the tRNA library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/tRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_tRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/tRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - 2> /dev/null | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_tRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to tRNAs of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## map ALL READS to piRNAs
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_piRNA.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to piRNAs of the primary organism:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that map to the genome to the piRNA library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/piRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_piRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/piRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_piRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that DO NOT map to the genome to the piRNA library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/piRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_piRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/piRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_piRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to piRNAs of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## map ALL READS to gencode transcripts
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_gencode.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to all ensembl/gencode transcripts of the primary organism:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that map to the genome to the gencode library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/gencode - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_gencode.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/gencode - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_gencode.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that DO NOT map to the genome to the gencode library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/gencode - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_gencode.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/gencode - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_gencode.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to all ensembl/gencode transcripts of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## map ALL READS to circular RNA junctions
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_circRNA.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to circular RNA junctions of the primary organism:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that map to the genome to the circular RNAlibrary
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/circularRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_circRNA.log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/circularRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_circRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Map reads that DO NOT map to the genome to the circular RNA library
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/circularRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_circRNA.log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/circularRNAs - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_circRNA.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to circular RNA junctions of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## process alignments
$(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped_noLibs.fq: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_miRNA_precursor.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_tRNA.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_piRNA.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_circRNA.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_gencode.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Collecting alignments...\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#$(SAMTOOLS_EXE) view $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam | awk -F "\t" '{print $$1"\t"$$2"\tunannotated:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_miRNA_precursor.bam | awk -F "\t" '{print $$1"\t"$$2"\tgenome:miRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_miRNA_precursor.bam | awk -F "\t" '{print $$1"\t"$$2"\tnogenome:miRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_tRNA.bam | awk -F "\t" '{print $$1"\t"$$2"\tgenome:tRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_tRNA.bam | awk -F "\t" '{print $$1"\t"$$2"\tnogenome:tRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_piRNA.bam | awk -F "\t" '{print $$1"\t"$$2"\tgenome:piRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_piRNA.bam | awk -F "\t" '{print $$1"\t"$$2"\tnogenome:piRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_gencode.bam | awk -F "\t" '{print $$1"\t"$$2"\tgenome:gencode:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14"\t"$$15}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_gencode.bam | awk -F "\t" '{print $$1"\t"$$2"\tnogenome:gencode:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14"\t"$$15}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeMapped_circRNA.bam | awk -F "\t" '{print $$1"\t"$$2"\tgenome:circRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_genomeUnmapped_circRNA.bam | awk -F "\t" '{print $$1"\t"$$2"\tnogenome:circRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished collecting alignments\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Sort alignments
	@echo -e "$(ts) $(SMRNAPIPELINE): Sorting alignments...\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam | sort -k 1,1 | sed 's/ /:/g' > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_LIBS.sam
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished sorting alignments\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Assign reads
	@echo -e "$(ts) $(SMRNAPIPELINE): Assigning reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) ProcessEndogenousAlignments --libPriority $(ENDOGENOUS_LIB_PRIORITY) --hairpin2genome $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRBase_v21_hairpin_hsa_hg19_aligned.sam --mature2hairpin $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRBase_v21_mature_hairpin_hsa_aligned.sam --reads2all $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_LIBS.sam --outputPath $(OUTPUT_DIR)/$(SAMPLE_ID) $(ENDOGENOUS_QUANT_RANDOM_BARCODE_STATS)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) ProcessEndogenousAlignments --libPriority $(ENDOGENOUS_LIB_PRIORITY) --hairpin2genome $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRNA_precursor2genome.sam --mature2hairpin $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/miRNA_mature2precursor.sam --reads2all $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_LIBS.sam --outputPath $(OUTPUT_DIR)/$(SAMPLE_ID) $(ENDOGENOUS_QUANT_RANDOM_BARCODE_STATS) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished assigning reads\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Summarise alignment statistics
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_miRNAmature_sense.txt | awk '{SUM+=$$4}END{printf "miRNA_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_miRNAmature_antisense.txt | awk '{SUM+=$$4}END{printf "miRNA_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_miRNAprecursor_sense.txt | awk '{SUM+=$$4}END{printf "miRNAprecursor_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_miRNAprecursor_antisense.txt | awk '{SUM+=$$4}END{printf "miRNAprecursor_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_tRNA_sense.txt | awk '{SUM+=$$4}END{printf "tRNA_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_tRNA_antisense.txt | awk '{SUM+=$$4}END{printf "tRNA_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_piRNA_sense.txt | awk '{SUM+=$$4}END{printf "piRNA_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_piRNA_antisense.txt | awk '{SUM+=$$4}END{printf "piRNA_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_gencode_sense.txt | awk '{SUM+=$$4}END{printf "gencode_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_gencode_antisense.txt | awk '{SUM+=$$4}END{printf "gencode_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_circularRNA_sense.txt | awk '{SUM+=$$4}END{printf "circularRNA_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_circularRNA_antisense.txt | awk '{SUM+=$$4}END{printf "circularRNA_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Count reads not mapping to the genome or to the libraries
	@echo -e "$(ts) $(SMRNAPIPELINE): Outputting reads not aligned to either the genome or transcripts:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	grep "nogenome" $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_Accepted.txt | awk '{print $$1}' | uniq > $(OUTPUT_DIR)/$(SAMPLE_ID)/readsMappedToLibs.tmp
	gzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_Accepted.txt > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_Accepted.txt.gz
	$(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) FilterFastxByIDList -b -p -IDs $(OUTPUT_DIR)/$(SAMPLE_ID)/readsMappedToLibs.tmp $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped.fq > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped_noLibs.fq 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	wc -l $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped_noLibs.fq | awk '{print "not_mapped_to_genome_or_libs\t"($$1)/4}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished outputting reads not aligned to either the genome or transcripts\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Tidy up
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/*.tmp
	#rm $(OUTPUT_DIR)/$(SAMPLE_ID)/readCounts_null.txt
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_Accepted.txt
	gzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_LIBS.sam > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_LIBS.sam.gz
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_LIBS.sam



##
## Align reads to repetitive element sequences, just in case repetitive reads have not been mapped to the genome
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotRepetitive.fq: $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped_noLibs.fq $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_unspliced.bam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to repetitive elements in the host genome:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped_noLibs.fq | $(BOWTIE2_EXE) -p $(N_THREADS) --local -D 20 -R 3 -N $(BOWTIE2_MAX_MISMATCHES) -L $(BOWTIE_SEED_LENGTH) -i S,1,0.50 -k $(N_GENOME_MULTIMAPS_TO_REPORT) --un $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotEndogenous.fq -x $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/repetitiveElements -U - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b -  > $(OUTPUT_DIR)/$(SAMPLE_ID)/RepeatElementsMapped.bam\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousUnaligned_ungapped_noLibs.fq | $(BOWTIE2_EXE) -p $(N_THREADS) --local -D 20 -R 3 -N $(BOWTIE2_MAX_MISMATCHES) -L $(BOWTIE_SEED_LENGTH) -i S,1,0.50 -k $(N_GENOME_MULTIMAPS_TO_REPORT) --un $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotRepetitive.fq -x $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/repetitiveElements -U - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/RepeatElementsMapped.bam
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to repetitive elements in the host genome\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Assigned non-redundantly to annotated REs	
	$(SAMTOOLS_EXE) view $(OUTPUT_DIR)/$(SAMPLE_ID)/RepeatElementsMapped.bam | grep -v "^@" | awk '{print $$1}' | sort | uniq | wc -l | awk '{print "repetitiveElements\t"$$0}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats



##
## map REMAINING reads to the genome (bowtie 2 gapped)
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotEndogenous.fq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotRepetitive.fq
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Aligning remaining reads to the genome allowing gaps \n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotRepetitive.fq | $(BOWTIE2_EXE) -p $(N_THREADS) --local -D 20 -R 3 -N $(BOWTIE2_MAX_MISMATCHES) -L $(BOWTIE_SEED_LENGTH) -i S,1,0.50 -k $(N_GENOME_MULTIMAPS_TO_REPORT) --reorder --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotEndogenous.fq.gz -x $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/$(MAIN_ORGANISM_GENOME_ID) -U - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_gapped.bam\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotRepetitive.fq | $(BOWTIE2_EXE) -p $(N_THREADS) --local -D 20 -R 3 -N $(BOWTIE2_MAX_MISMATCHES) -L $(BOWTIE_SEED_LENGTH) -i S,1,0.50 -k $(N_GENOME_MULTIMAPS_TO_REPORT) --reorder --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotEndogenous.fq.gz -x $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/$(MAIN_ORGANISM_GENOME_ID) -U - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_gapped.bam
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished aligning remaining reads to the genome allowing gaps\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## mapped to the genome with gaps
	$(SAMTOOLS_EXE) view $(OUTPUT_DIR)/$(SAMPLE_ID)/endogenousAlignments_gapped.bam | grep -v "^@" | awk '{print $$1}' | sort | uniq | wc -l | awk '{print "endogenous_gapped\t"$$0}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats



##
## Use the unmapped reads and search against all miRNAs in miRBase
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/unaligned.fq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotEndogenous.fq.gz
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to all miRNAs in miRBase:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotEndogenous.fq.gz | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata --un $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/unaligned.fq $(DATABASE_PATH)/miRBase/miRBase_precursors - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.bowtie1stats | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/reads_NotEndogenous.fq.gz | $(BOWTIE1_EXE) -p $(N_THREADS) --chunkmbs $(BOWTIE_CHUNKMBS) -l $(BOWTIE_SEED_LENGTH) -n $(BOWTIE1_MAX_MISMATCHES) --all --sam --fullref --best --strata --un $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/unaligned.fq $(DATABASE_PATH)/miRBase/miRBase_precursors - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.bowtie1stats | awk '{if($$1 ~ /^@/ || $$2 != 4){print $$0}}' | $(SAMTOOLS_EXE) view -@ $(N_THREADS) -b - | $(SAMTOOLS_EXE) sort -n -m $(SAMTOOLS_SORT_MEM) -@ $(N_THREADS) -O bam -T $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/tmp - > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.bam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	gzip $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/unaligned.fq
	## quantify read alignments using a slight hack of the endogenous alignment engine
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.bam | sort -k1,1 | awk -F "\t" '{print $$1"\t"$$2"\tnogenome:miRNA:"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"$$9"\t"$$10"\t"$$11"\t"$$12"\t"$$13"\t"$$14}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.sam
	@echo -e "$(ts) $(SMRNAPIPELINE): Assigning reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) ProcessEndogenousAlignments --libPriority miRNA --hairpin2genome $(DATABASE_PATH)/miRBase/miRNA_precursor2genome.sam --mature2hairpin $(DATABASE_PATH)/miRBase/miRNA_mature2precursor.sam --reads2all $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.sam --outputPath $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) ProcessEndogenousAlignments --libPriority miRNA --hairpin2genome $(DATABASE_PATH)/miRBase/miRNA_precursor2genome.sam --mature2hairpin $(DATABASE_PATH)/miRBase/miRNA_mature2precursor.sam --reads2all $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.sam --outputPath $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished assigning reads\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Tidy up:
	gzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/endogenousAlignments_Accepted.txt > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenousMiRNAAlignments_Accepted.txt.gz
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/endogenousAlignments_Accepted.txt
	## Input to exogenous miRNA alignment
	grep "# reads processed:" $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.bowtie1stats | awk -F ": " '{print "input_to_exogenous_miRNA\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Assigned non-redundantly to annotated exogenous miRNAs
	grep "# reads with at least one reported alignment:" $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/exogenous_miRBase_mapped.bowtie1stats | awk -F ": " '{print $$2}' | awk -F " " '{print "exogenous_miRNA\t"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to all miRNAs in miRBase\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log



##
## Use the unmapped reads and search against all rRNAs in RDP (ribosome DB)
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/unaligned.fq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/unaligned.fq.gz
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to rRNA sequences in RDP (a.k.a. ribosome DB):\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): $(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_ --genomeDir $(DATABASE_PATH)/ribosomeDatabase/exogenous_rRNAs --readFilesIn $< --outReadsUnmapped Fastx --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_ --genomeDir $(DATABASE_PATH)/ribosomeDatabase/exogenous_rRNAs --readFilesIn $< --outReadsUnmapped Fastx --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Input to exogenous rRNA alignment
	grep "Number of input reads" $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Log.final.out | tr '[:blank:]' ' ' | awk -F " \\\| " '{print "input_to_exogenous_rRNA\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Assigned non-redundantly to annotated exogenous rRNAs
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Aligned.out.sam | grep -v "^@" | awk '{print $$1}' | sort | uniq | wc -l | awk '{print "exogenous_rRNA\t"$$0}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## compress and tidy up
	gzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Unmapped.out.mate1 > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/unaligned.fq.gz
	$(SAMTOOLS_EXE) view -@ $(N_THREADS) -b $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Aligned.out.sam > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Aligned.bam
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Unmapped.out.mate1
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Aligned.out.sam
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/exogenous_rRNA_Log.out
	@echo -e "$(ts) $(SMRNAPIPELINE): Finished mapping to rRNA sequences in RDP\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log



##
## Routines for aligning unmapped reads to exogenous sequences
##
## Bacteria
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria10_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/unaligned.fq.gz
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to exogenous GENOMES of BACTERIA:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Bacteria1:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria1_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA1 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria2:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria2_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA2 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria3:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria3_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA3 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria4:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria4_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA4 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria5:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria5_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA5 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria6:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria6_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA6 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria7:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria7_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA7 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria8:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria8_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA8 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria9:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria9_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA9 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Bacteria10:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria10_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA10 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Finished mapping to exogenous GENOMES of BACTERIA\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## Plants
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/unaligned.fq.gz
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to exogenous GENOMES of PLANTS:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Plants1:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants1_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS1 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Plants2:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants2_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS2 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Plants3:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants3_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS3 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Plants4:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants4_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS4 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Plants5:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS5 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Finished mapping to exogenous GENOMES of PLANTS\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## Fungi, Protist, and Virus
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/unaligned.fq.gz
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to exogenous GENOMES of FUNGI, PROTISTS, and VIRUSES:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_FUNGI_PROTIST_VIRUS --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Finished mapping to exogenous GENOMES of FUNGI, PROTISTS, and VIRUSES\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

## Vertebrates
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/unaligned.fq.gz
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Mapping reads to exogenous GENOMES of VERTEBRATES:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Vertebrate1:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate1_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE1 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Vertebrate2:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate2_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE2 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Vertebrate3:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate3_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE3 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Vertebrate4:" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE4 --readFilesIn $< --parametersFiles $(STAR_PARAMS_FILE_PATH) --outFilterMismatchNoverLmax $(STAR_outFilterMismatchNoverLmax) --outFilterMismatchNmax $(BOWTIE2_MAX_MISMATCHES) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Finished mapping to exogenous GENOMES of VERTEBRATES\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Combine exogenous genome alignment info
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria10_Aligned.out.sam $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) $(SMRNAPIPELINE): Collecting all alignments to exogenous GENOMES:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#
	#cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria1_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | sed 's/[:$$]/ /g' | awk '{print $$1"\tBacteria\t"$$2"\t"$$7"\t"$$12"\t"$$13"\t"$$14}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria1_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria2_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria3_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria4_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria5_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria6_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria7_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria8_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria9_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria10_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[$$:]");print $$1"\tBacteria\t"a[1]"\t"a[7]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	#
	#cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants1_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants1_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants2_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants3_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants4_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants_Aligned.out.sam.summary
	#
	#cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam | grep -v "^@" | grep "Virus:" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | sed 's/|/ /g' | awk '{print $$1"\t"$$2"\t"$$4"\t"$$6"\t"$$7"\t"$$8"\t"$$9}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Virus_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam | grep -v "^@" | grep "Virus:" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,"[:|]");print $$1"\t"a[1]"\t"a[3]"\t"a[5]"\t"$$3"\t"$$4"\t"$$5}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Virus_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam | grep -v "^@" | grep "Fungi:" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Fungi_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam | grep -v "^@" | grep "Protist:" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Protist_Aligned.out.sam.summary
	#
	#cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate1_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate1_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate2_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate3_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam | grep -v "^@" | awk '{print $$1,$$3,$$4,$$6,$$10}' | uniq | awk '{split($$2,a,":");print $$1"\t"a[1]"\t"a[2]"\t"a[3]"\t"$$3"\t"$$4"\t"$$5}' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate_Aligned.out.sam.summary
	#
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Virus_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Fungi_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Protist_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	#
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt | sort -k 1,1 > $@
	#
	## Input to exogenous genome alignment
	gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_rRNA/unaligned.fq.gz | wc -l | awk '{print "input_to_exogenous_genomes\t"$$1/4}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Count reads mapped to exogenous genomes:
	#cat $@ | awk '{print $$1}' | uniq | awk -F "#" '{SUM += $$2} END {print "exogenous_genomes\t"SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $@ | awk '{print $$1}' | uniq | wc -l | awk '{print "exogenous_genomes\t"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	#
	@echo -e "\n$(ts) $(SMRNAPIPELINE): Finished mapping to exogenous GENOMES of VERTEBRATES\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Create exogenous alignment result matrix:
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.result.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.txt
	# remove duplicate reads that correspond to multimaps within a single species
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.txt | awk '{print $$1"\t"$$2"\t"$$3}' | uniq > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.unique.txt
	# get the IDs of reads aligning uniquely to one kingdom
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.unique.txt | awk '{print $$1"\t"$$2}' | uniq | awk '{print $$1}' | uniq -c | awk '{if($$1==1) print $$2}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.kingdomUniqueReads
	# get the IDs of reads aligning uniquely to one species
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.unique.txt | awk '{print $$1}' | uniq -c | awk '{if($$1==1) print $$2}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.speciesUniqueReads
	#
	## count all reads aligning to each species, regardless of multimapping
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.unique.txt | awk '{arr[$$2"\t"$$3]+=1} END {for (x in arr) print x"\t"arr[x]}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.counts
	## count only reads aligning to each species that only multi map to the same kingdom
	awk 'NR==FNR {h[$$1]="YES_YES_YES"; next} {print $$1,$$2,$$3,h[$$1]}' $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.kingdomUniqueReads $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.unique.txt | grep "YES_YES_YES" | awk '{print $$1"\t"$$2"\t"$$3}' | awk '{arr[$$3]+=1} END {for (x in arr) print x"\t"arr[x]}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.kingdom.counts
	awk 'NR==FNR {h[$$1]="YES_YES_YES"; next} {print $$1,$$2,$$3,h[$$1]}' $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.speciesUniqueReads $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.unique.txt | grep "YES_YES_YES" | awk '{print $$1"\t"$$2"\t"$$3}' | awk '{arr[$$3]+=1} END {for (x in arr) print x"\t"arr[x]}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.species.counts
	## combine all counts and kingdom level counts:
	awk 'NR==FNR {h[$$1]=$$2; next} {if($$2 in h) print $$0"\t"h[$$2]; else print $$0"\t0"}' $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.kingdom.counts $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.counts > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp
	#
	## output results
	@echo -e "Kingdom\tSpecies\tReadCount_allMappedToThisSpecies\tReadCount_kingdomSpecific\tReadCount_speciesSpecific" > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.result.txt
	awk 'NR==FNR {h[$$1]=$$2; next} {if($$2 in h) print $$0"\t"h[$$2]; else print $$0"\t0"}' $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp.species.counts $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp | sort -nrk 5,5 >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.result.txt
	#
	# Final output format of ExogenousGenomicAlignments.result.txt:
	# Kingdom	Species	ReadCount_allMappedToThisSpecies	ReadCount_kingdomSpecific	ReadCount_speciesSpecific
	# Plant	Oryza_indica	25413	2334	559
	#
	# tidy up
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/tmp*




##
##
##