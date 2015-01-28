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
## Version 2.0.6 (2015-01-27)                                                        ##
##                                                                                   ##
#######################################################################################


##
## 1) On the command line, be sure to specify the following MANDATORY parameters
##
DATA_DIR                := NULL
OUTPUT_DIR              := NULL
INPUT_FILE_PATH         := NULL
SAMPLE_NAME             := NULL
## You can also override the following OPTIONAL parameters on the commandline
CALIBRATOR_LIBRARY      := NULL


##
## 2) Choose the main organism for smallRNA / genome alignment (hsa + hg19, hsa + hg38, or mmu + mm10)
##
## Human:
MAIN_ORGANISM           := hsa
MAIN_ORGANISM_GENOME_ID := hg38
## Mouse:
#MAIN_ORGANISM           := mmu
#MAIN_ORGANISM_GENOME_ID := mm10


##
## 3) Select whether pipeline is run locally, should be 'true' unless this is the Genboree implementation!
##
LOCAL_EXECUTION := true


##
## 4) Choose optional analysis-specific options (or specify these at the command line)
##
ADAPTER_SEQ                     := NULL
TRNA_MAPPING                    := on
PIRNA_MAPPING                   := on
GENCODE_MAPPING                 := on
REPETITIVE_ELEMENT_MAPPING      := on
REMOVE_LARGE_INTERMEDIATE_FILES := false

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
	N_THREADS := 4
	JAVA_RAM  := 12G
	MAX_RAM   := 12000000000
	## NB: The 'EXE_DIR' MUST be an ABSOLUTE PATH or sRNABench will fail!
	EXE_DIR   := /gpfs/scratch/fas/gerstein/rrk24/bin/smallRNAPipeline
	
	
	##
	## 6) Check that the paths to the required 3rd party executables work!
	##
	JAVA_EXE         := /usr/bin/java
	#FASTX_ARTFCT_EXE := $(EXE_DIR)/fastx_0.0.14/bin/fastx_artifacts_filter
	FASTX_CLIP_EXE   := $(EXE_DIR)/fastx_0.0.14/bin/fastx_clipper
	FASTX_FILTER_EXE := $(EXE_DIR)/fastx_0.0.14/bin/fastq_quality_filter
	#BOWTIE1_PATH    := $(EXE_DIR)/bowtie-1.1.1
	#VIENNA_PATH     := $(EXE_DIR)/ViennaRNA_2.1.5/bin
	BOWTIE_EXE       := $(EXE_DIR)/bowtie2-2.2.4/bowtie2
	SAMTOOLS_EXE     := $(EXE_DIR)/samtools-0.1.18/samtools
	FASTQC_EXE       := $(JAVA_EXE) -classpath $(EXE_DIR)/FastQC_0.11.2:$(EXE_DIR)/FastQC_0.11.2/sam-1.103.jar:$(EXE_DIR)/FastQC_0.11.2/jbzip2-0.9.jar
	SRATOOLS_EXE     := $(EXE_DIR)/sratoolkit.2.1.7-centos_linux64/fastq-dump
	THUNDER_EXE      := $(EXE_DIR)/Thunder.jar
	SRNABENCH_EXE    := $(EXE_DIR)/sRNAbench.jar
	SRNABENCH_LIBS   := $(EXE_DIR)/sRNAbenchDB
	STAR_EXE         := $(EXE_DIR)/STAR_2.4.0i/bin/Linux_x86_64/STAR
	STAR_GENOMES_DIR := /gpfs/scratch/fas/gerstein/rrk24/ANNOTATIONS/Genomes_BacteriaFungiMammalPlantProtistVirus
	
	##
	## Use the input path to infer filetype and short name
	##
	INPUT_FILE_NAME := $(notdir $(INPUT_FILE_PATH))
	INPUT_FILE_ID   := $(basename $(INPUT_FILE_NAME))
	
else
	##
	## These parameters are for the Genboree installation only
	##
	EXE_DIR := $(SCRATCH_DIR)
	N_THREADS := $(N_THREADS)
	JAVA_RAM := 12G
	MAX_RAM := 12000000000

	#FASTX_ARTFCT_EXE := fastx_artifacts_filter
	FASTX_CLIP_EXE := fastx_clipper
	FASTX_FILTER_EXE := fastq_quality_filter
	BOWTIE1_PATH := NULL
	VIENNA_PATH := NULL
	BOWTIE_EXE := bowtie2
	SAMTOOLS_EXE := samtools
	FASTQC_EXE := $(JAVA_EXE) -classpath fastqc
	SRATOOLS_EXE := fastq-dump
	SRNABENCH_EXE := $(SRNABENCH_EXE)
	THUNDER_EXE := $(THUNDER_EXE)
	
	## Path to sRNABench libraries
	SRNABENCH_LIBS := $(SRNABENCH_LIBS)
	
	INPUT_FILE_NAME := $(notdir $(INPUT_FILE_PATH))
    INPUT_FILE_ID := $(INPUT_FILE_ID)
endif



## Define current time
ts := `/bin/date "+%Y-%m-%d--%H:%M:%S"`

##
## Initialise organism specific smallRNA library parameters
##
## override format: #<N_mismatches>#<seed_length>#<alignMode [v|n]>#<N_multimaps>
MISMATCH_N_MIRNA := 1
MISMATCH_N_OTHER := 2
MULTIMAP_MAX     := 10000
BOWTIE_OVERRIDE  := \#$(MISMATCH_N_OTHER)\#9\#v\#$(MULTIMAP_MAX)

ifeq ($(MAIN_ORGANISM),hsa)  ## FOR HUMAN
	
	## Override the genome for adapter identification (saves having bt2 indexes for both hg19 and hg38)
	GENOME_ID_FOR_ADAPTER := hg19
	INDEX_TRNA			  := hg19_tRNAs

	ifeq ($(MAIN_ORGANISM_GENOME_ID),hg19) ## hg19

		INDEX_GENCODE 		:= hg19_gencode
		INDEX_PIRNA	 		:= hg19_piRNAs
		INDEX_REP           := hg19_RepetitiveElements
		BOWTIE_INDEX_RRNA	:= $(SRNABENCH_LIBS)/customIndices/hg19_rRNA
		
	else ifeq ($(MAIN_ORGANISM_GENOME_ID),hg38) ## hg38
		
		INDEX_GENCODE 		:= hg38_gencode
		INDEX_PIRNA	 		:= hg19_piRNAs
		INDEX_REP           := hg38_RepetitiveElements
		BOWTIE_INDEX_RRNA	:= $(SRNABENCH_LIBS)/customIndices/hg38_rRNA
		
	endif

else ifeq ($(MAIN_ORGANISM),mmu)  ## FOR MOUSE
	
	GENOME_ID_FOR_ADAPTER := mm10

	INDEX_GENCODE 		:= mm10_gencode
	INDEX_TRNA			:= mm10_tRNAs
	INDEX_PIRNA	 		:= mm10_piRNAs
	INDEX_REP           := mm10_RepetitiveElements
	BOWTIE_INDEX_RRNA	:= $(SRNABENCH_LIBS)/customIndices/mm10_rRNA
	
endif

GENCODE_LIBS := libs=$(INDEX_GENCODE)$(BOWTIE_OVERRIDE)
TRNA_LIBS    := libs=$(INDEX_TRNA)$(BOWTIE_OVERRIDE)
PIRNA_LIBS   := libs=$(INDEX_PIRNA)$(BOWTIE_OVERRIDE)
REP_LIBS     := libs=$(INDEX_REP)$(BOWTIE_OVERRIDE)


##
## Turn off a smallRNA library if they are not selected by the user
##
ifneq ($(TRNA_MAPPING),on)
	TRNA_LIBS    :=
	INDEX_TRNA   :=
endif
ifneq ($(PIRNA_MAPPING),on)
	PIRNA_LIBS   :=
endif
ifneq ($(GENCODE_MAPPING),on)
	GENCODE_LIBS :=
endif
ifneq ($(REPETITIVE_ELEMENT_MAPPING),on)
	REP_LIBS     :=
endif




USEAGE := 
ifeq ($(INPUT_FILE_ID),NULL)
  USEAGE := "make -f smallRNA_pipeline INPUT_FILE_PATH=[required: absolute/path/to/input/.fa|.fq|.sra] N_THREADS=[required: number of threads] OUTPUT_DIR=<required: absolute/path/to/output> INPUT_FILE_ID=[required: samplename] ADAPTER_SEQ=[optional: will guess sequence if not provided here; none, if already clipped input] MAIN_ORGANISM=[optional: defaults to 'hsa'] MAIN_ORGANISM_GENOME_ID=[optional: defaults to 'hg38'] CALIBRATOR_LIBRARY=[optional: path/to/bowtie/index/containing/calibrator/sequences] TRNA_MAPPING=[optional: TRUE|FALSE, default is TRUE] GENCODE_MAPPING=[optional: TRUE|FALSE, default is TRUE] PIRNA_MAPPING=[optional: TRUE|FALSE, default is TRUE] MAP_EXOGENOUS=[optional: off|miRNA|on, default is miRNA]"
endif



##
## Try to export the Bowtie1 and ViennaRNA executable directories to the PATH
##
EXPORT_CMD :=
ifneq ($(BOWTIE1_EXE),NULL)
	EXPORT_CMD := export PATH=$$PATH:$(BOWTIE1_PATH):$(VIENNA_PATH)
endif


## Path to genome bowtie2 index
BOWTIE_INDEX_GENOME := $(SRNABENCH_LIBS)/index/$(GENOME_ID_FOR_ADAPTER)


## Path to the UniVec contaminants DB
BOWTIE_INDEX_UNIVEC := $(SRNABENCH_LIBS)/customIndices/UniVec_Core.contaminants


## SmallRNA sequence libraries to map against AFTER mapping to the known miRNAs for the target organism (see below)
OTHER_LIBRARIES := $(TRNA_LIBS) $(PIRNA_LIBS) $(GENCODE_LIBS) $(REP_LIBS)


## Map reads to plant and virus miRNAs
ifeq ($(MAP_EXOGENOUS),miRNA)		## ALIGNMENT TO ONLY EXOGENOUS MIRNA
	PROCESS_SAMPLE_REQFILE := EXOGENOUS_miRNA/reads.fa
else ifeq ($(MAP_EXOGENOUS),on)	## COMPLETE EXOGENOUS GENOME ALIGNMENT
	PROCESS_SAMPLE_REQFILE := EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.txt
else
	PROCESS_SAMPLE_REQFILE := reads.fa
endif


##
## List of plant and virus species IDs to which to map reads that do not map to the genome of the primary organism
##
EXOGENOUS_MIRNA_SPECIES := $(shell cat $(SRNABENCH_LIBS)/libs/mature.fa | grep ">" | awk -F '-' '{print $$1}' | sed 's/>//g'| sort | uniq | tr '\n' ':' | rev | cut -c 2- | rev | sed 's/$(MAIN_ORGANISM)://g')

## Parameters to use for the bowtie mapping of calibrator oligos and rRNAs
BOWTIE2_MAPPING_PARAMS_CALIBRATOR := -D 15 -R 2 -N 1 -L 16 -i S,1,0
BOWTIE2_MAPPING_PARAMS_RRNA       := -D 15 -R 2 -N 1 -L 19 -i S,1,0

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
#$(COMMAND_CONVERT_SRA) | head -n 4000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "Phred_33"; else if(max>73 && min>=64) print "Phred_64"; else if(min>=59 && min<64 && max>73) print "Solexa_64"; else print "Unknown";}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding
COMMAND_GUESS_QUALITY_ENCODING := $(COMMAND_CONVERT_SRA) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding


##
## Logic block to write the adapter sequence (whether or not one is provided by the user) to the .adapterSeq file
##
ifeq ($(ADAPTER_SEQ),NULL)
	COMMAND_WRITE_ADAPTER_SEQ := $(COMMAND_CONVERT_SRA) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 40000000 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | $(BOWTIE_EXE) --no-head -p $(N_THREADS) --local -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 -k 2 --upto 10000000 -x $(BOWTIE_INDEX_GENOME) -U - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{if ($$5==255) print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.unique.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.unique.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{print $$6}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort -rnk 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 100 > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).cigarFreqs 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.unique.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{if ($$2==0) print $$3"\t"$$4"\t"$$6"\t"$$10}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | grep "[[:space:]]2[0-9]M[0-9][0-9]S" > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{print $$3}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort -rnk 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 100 > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).okCigarFreqs 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).okCigarFreqs 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | head -n 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{print substr($$2,1,2)}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.txt 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.sam 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | grep "[[:space:]]$$(<$(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.txt)" 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{getline len<"$(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.txt"; print substr($$4,len+1)}' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sed 's/[A]*$$//' 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | sort -rnk 1 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | awk '{if ($$1 > 75) print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).potentialAdapters.txt 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err;  \
	head -n 1 $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).potentialAdapters.txt | awk '{print $$2}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq;  \
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/tmp.*
	LOGENTRY_WRITE_ADAPTER := $(ts) SMRNAPIPELINE: Removing 3' adapter sequence using fastX:\n
else ifeq ($(ADAPTER_SEQ),none)
	COMMAND_WRITE_ADAPTER_SEQ := echo 'no adapter' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq;
	COMMAND_CLIP_ADAPTER := $(COMMAND_CONVERT_SRA) | gzip -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	LOGENTRY_WRITE_ADAPTER := Provided 3' adapter clipped input sequence file. No clipping necessary.\n 
else
	COMMAND_WRITE_ADAPTER_SEQ := echo $(ADAPTER_SEQ) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq
	LOGENTRY_WRITE_ADAPTER := $(ts) SMRNAPIPELINE: Provided 3' adapter sequence. Removing 3' adapter sequence using fastX:\n
endif


## If no adapter clipping command has been set- use this one:
COMMAND_CLIP_ADAPTER ?= $(COMMAND_CONVERT_SRA) | $(FASTX_CLIP_EXE) -a $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq) -l 10 -v -M 7 -z -o $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err




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
	COMMAND_MAP_CALIBRATOR := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_CALIBRATOR) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noCalibrator.fastq.gz -x $(CALIBRATOR_LIBRARY) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | tee $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.calibratormapped.bam | $(SAMTOOLS_EXE) view - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.calibratormapped.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	LOGENTRY_MAP_CALIBRATOR_1 := $(ts) SMRNAPIPELINE: Mapping reads to calibrator sequences using bowtie:\n
	LOGENTRY_MAP_CALIBRATOR_2 := $(ts) SMRNAPIPELINE: Finished mapping to the calibrators\n
	
	FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT := $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noCalibrator.fastq.gz
	
endif


##
## Bowtie2 command to align reads to the UniVec contaminant sequence database
##
COMMAND_MAP_UNIVEC := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noUniVecContaminants.fastq.gz -x $(BOWTIE_INDEX_UNIVEC) -U $(FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam | grep -v "^@" | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminants.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam | grep -v "^@" | awk '{print $$1}' | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c | wc -l > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminants.readCount 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.bam; \
rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminantMapped.sam


##
## Bowtie2 command to align reads to the rRNA sequences
##
COMMAND_MAP_RRNAS := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz -x $(BOWTIE_INDEX_RRNA) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noUniVecContaminants.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam | grep -v "^@" | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNA.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam | grep -v "^@" | awk '{print $$1}' | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c | wc -l > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNA.readCount 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.bam; \
rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sam



## Parameters for the endogenous-exRNA mapping
FIXED_PARAMS_MAIN      := -Xmx$(JAVA_RAM) -jar $(SRNABENCH_EXE) dbPath=$(SRNABENCH_LIBS) p=$(N_THREADS) chunkmbs=2000 microRNA=$(MAIN_ORGANISM) species=$(MAIN_ORGANISM_GENOME_ID) plotMiR=true plotLibs=false predict=false $(OTHER_LIBRARIES) writeGenomeDist=true noMM=$(MISMATCH_N_MIRNA) maxReadLength=75 noGenome=true tRNA=$(INDEX_TRNA) mBowtie=$(MULTIMAP_MAX)
## Parameters for the exogenous-exRNA mapping
FIXED_PARAMS_EXOGENOUS := -Xmx$(JAVA_RAM) -jar $(SRNABENCH_EXE) dbPath=$(SRNABENCH_LIBS) p=$(N_THREADS) chunkmbs=2000 microRNA=$(EXOGENOUS_MIRNA_SPECIES) plotMiR=true predict=false noMM=$(MISMATCH_N_MIRNA)


##
## Remove some potentially large intermediate pipeline output (can save as much as 50% total output size)
##
TIDYUP_COMMAND := 
ifeq ($(REMOVE_LARGE_INTERMEDIATE_FILES),true)
	TIDYUP_COMMAND := rm $(OUTPUT_DIR)/$(SAMPLE_ID)/genome.parsed; rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped*.fastq.gz
endif




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
## Delete sample results and logfiles
##
#clean: 
#	rm -r $(OUTPUT_DIR)/$(SAMPLE_ID)




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
	@echo -e "$(ts) SMRNAPIPELINE: BEGIN smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" > $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: BEGIN \n" > $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Processing adapter sequence:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_WRITE_ADAPTER_SEQ)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_WRITE_ADAPTER_SEQ)
	@echo -e "$(ts) SMRNAPIPELINE: Progress_1_FoundAdapter" > $(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat
	#
	@echo -e "#STATS from smallRNA-seq Pipeline for sample $(SAMPLE_ID)" > $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	@echo -e "Stage\tReadCount" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## CLIP 3' adapter sequence
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Adapter sequence: $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(LOGENTRY_WRITE_ADAPTER)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_CLIP_ADAPTER)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_CLIP_ADAPTER)
	@echo -e "$(ts) SMRNAPIPELINE: Finished removing adapters\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count reads input to adapter clipping
	grep "Input: " $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "input\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Count reads output following adapter clipping
	grep "Output: " $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "clipped\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## FILTER clipped reads that have poor overall base quality
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Guessing encoding of fastq read-qualities:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_GUESS_QUALITY_ENCODING)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_GUESS_QUALITY_ENCODING)
	@echo -e "$(ts) SMRNAPIPELINE: Finished guessing encoding of fastq read-qualities:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Filtering reads by base quality:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: gunzip -c $< | $(FASTX_FILTER_EXE) -vz -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $@\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	gunzip -c $< | $(FASTX_FILTER_EXE) -vz -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $@ 2>>$(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Finished filtering reads by base quality\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count reads that failed the quality filter
	grep "low-quality reads" $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "failed_quality_filter\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

##
## Assess Read-lengths after clipping
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.readLengths.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz
	@echo -e "======================" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Calculating length distribution of clipped reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) GetSequenceLengths $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq > $@ 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	gunzip -c $< > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq
	$(JAVA_EXE) -Xmx$(JAVA_RAM) -jar $(THUNDER_EXE) GetSequenceLengths $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq > $@ 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq
	@echo -e "$(ts) SMRNAPIPELINE: Finished calculating read-lengths\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Perform FastQC after adapter removal
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered_fastqc.zip: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz
	@echo -e "======================" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Running FastQC on clipped reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/bbsrc/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: Finished running FastQC\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## MAP to external bowtie (calibrator?) library and to the rRNA sequences
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.fastq.gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.readLengths.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered_fastqc.zip
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(LOGENTRY_MAP_CALIBRATOR_1)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_MAP_CALIBRATOR)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_CALIBRATOR)
	@echo -e "$(ts) SMRNAPIPELINE: $(LOGENTRY_MAP_CALIBRATOR_2)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count calibrator oligo reads
	$(COMMAND_COUNT_CALIBRATOR)
	#
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to contaminant sequences in UniVec using Bowtie2:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_MAP_UNIVEC)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_UNIVEC)
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to the UniVec contaminant DB\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count UniVec contaminant reads
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.uniVecContaminants.readCount | awk '{print "UniVec_contaminants\t"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	#
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to ribosomal RNA sequences using Bowtie2:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_MAP_RRNAS)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_RRNAS) 
	$(SAMTOOLS_EXE) sort $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sorted
	$(SAMTOOLS_EXE) index $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.sorted.bam
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNAmapped.bam
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to the rRNAs\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count rRNA reads
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.rRNA.readCount | awk ' {print "rRNA\t"$$1}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## Perform FastQC again after rRNA / UniVec removal
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA_fastqc.zip: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz
	@echo -e "======================" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Running FastQC on cleaned reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/bbsrc/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(FASTQC_EXE) -Xmx$(JAVA_RAM) -Dfastqc.threads=$(N_THREADS) -Dfastqc.unzip=false -Dfastqc.output_dir=$(OUTPUT_DIR)/$(SAMPLE_ID)/ uk/ac/babraham/FastQC/FastQCApplication $< >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: Finished running FastQC\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Map reads to the main genome of interest
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/noGenome/reads.fa: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA.fastq.gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.filtered.noRiboRNA_fastqc.zip
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to smallRNAs of primary organism:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(JAVA_EXE) $(FIXED_PARAMS_MAIN) input=$< output=$(OUTPUT_DIR)/$(SAMPLE_ID) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	$(JAVA_EXE) $(FIXED_PARAMS_MAIN) input=$< output=$(OUTPUT_DIR)/$(SAMPLE_ID) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to the small-RNAs of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: grep "chrUn_gl000220" $(OUTPUT_DIR)/$(SAMPLE_ID)/genome.txt | awk '{print $$1}' | sed 's/#/ /g' | awk '{ sum += $$2 } END { print "Number of reads mapped to chrUn_gl000220 = "sum }' >> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	grep "chrUn_gl000220" $(OUTPUT_DIR)/$(SAMPLE_ID)/genome.txt | awk '{print $$1}' | sed 's/#/ /g' | awk '{ sum += $$2 } END { print "Number of reads mapped to chrUn_gl000220 = "sum }' >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Count reads not mapped to rRNA
	grep "No. raw input reads:" $(OUTPUT_DIR)/$(SAMPLE_ID).log | head -n 1 | awk -F ':' '{print "not_rRNA\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Count reads mapped to the genome
	grep "out of" $(OUTPUT_DIR)/$(SAMPLE_ID)/summary.txt | awk '{print "genome\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Assigned non-redundantly to annotated miRNAs
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/mature_sense.grouped | awk '{sum+=$$4} END {printf "miRNA_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/mature_antisense.grouped | awk '{sum+=$$4} END {printf "miRNA_antisense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Assigned non-redundantly to annotated tRNAs
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_TRNA)_sense.grouped | awk '{sum+=$$4} END {printf "tRNA_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_TRNA)_antisense.grouped | awk '{sum+=$$4} END {printf "tRNA_antisense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Assigned non-redundantly to annotated piRNAs
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_PIRNA)_sense.grouped | awk '{sum+=$$4} END {printf "piRNA_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_PIRNA)_antisense.grouped | awk '{sum+=$$4} END {printf "piRNA_antisense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Assigned non-redundantly to annotated transcripts in Gencode
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_GENCODE)_sense.grouped | awk '{sum+=$$4} END {printf "Gencode_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_GENCODE)_antisense.grouped | awk '{sum+=$$4} END {printf "Gencode_antisense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	## Assigned non-redundantly to annotated repetitive elements
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_REP)_sense.grouped | awk '{sum+=$$4} END {printf "repetitiveElement_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_REP)_antisense.grouped | awk '{sum+=$$4} END {printf "repetitiveElement_antisense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## Use the unmapped reads and search against all plant and viral miRNAs in miRBase
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa: $(OUTPUT_DIR)/$(SAMPLE_ID)/noGenome/reads.fa
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to smallRNAs of plants and viruses:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(JAVA_EXE) $(FIXED_PARAMS_EXOGENOUS) input=$(OUTPUT_DIR)/$(SAMPLE_ID)/readsNotDetected.fa output=$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(JAVA_EXE) $(FIXED_PARAMS_EXOGENOUS) input=$(OUTPUT_DIR)/$(SAMPLE_ID)/noGenome/reads.fa output=$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to plant and virus small-RNAs\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	## Assigned non-redundantly to annotated exogenous miRNAs
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/mature_sense.grouped | awk '{sum+=$$4} END {printf "miRNA_exogenous_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## NEW routines for aligning unmapped reads to exogenous sequences
##
## Bacteria
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_BACTERIA --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in

## Plants
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants1_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS1 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants2_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS2 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants3_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS3 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants4_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS4 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_PLANTS5 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in

## Fungi, Protist, and Virus
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_FUNGI_PROTIST_VIRUS --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in

## Vertebrates
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate1_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE1 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate2_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE2 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate3_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE3 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in
	$(STAR_EXE) --runThreadN $(N_THREADS) --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_ --genomeDir $(STAR_GENOMES_DIR)/STAR_GENOME_VERTEBRATE4 --readFilesIn $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_miRNA/reads.fa --parametersFiles $(STAR_GENOMES_DIR)/STAR_Parameters_Exogenous.in



##
## Combine exogenous genome alignment info
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.sorted.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | sed 's/|/ /g' | awk '{print $$1"\t"$$2"\t"$$4"\t"$$6"\t"$$7"\t"$$8"\t"$$9}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary
	#
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants1_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants1_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants2_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants2_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants3_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants3_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants4_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants4_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam.summary
	#
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam | grep -v "^@" | grep "Virus:" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | sed 's/|/ /g' | awk '{print $$1"\t"$$2"\t"$$4"\t"$$6"\t"$$7"\t"$$8"\t"$$9}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Virus_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam | grep -v "^@" | grep "Fungi:" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Fungi_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/FungiProtistVirus_Aligned.out.sam | grep -v "^@" | grep "Protist:" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Protist_Aligned.out.sam.summary
	#
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate1_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate1_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate2_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate2_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate3_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate3_Aligned.out.sam.summary
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam | grep -v "^@" | awk '{print $$1" "$$3" "$$4" "$$6" "$$10}' | uniq | sed 's/:/ /g' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam.summary
	#
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Bacteria_Aligned.out.sam.summary > $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants1_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants2_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants3_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants4_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Plants5_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Virus_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Fungi_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Protist_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate1_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate2_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate3_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/Vertebrate4_Aligned.out.sam.summary >> $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt
	#
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/EXOGENOUS_genomes/ExogenousGenomicAlignments.txt | sort -k 1 > $@
	#
	## Count reads mapped to exogenous genomes:
	cat $@ | awk '{print $$1}' | uniq | awk -F "#" '{SUM += $$2} END {print "exogenous_genomes\t"SUM}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## Main sub-target
##
processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(PROCESS_SAMPLE_REQFILE)
	## Copy Output descriptions file
	cp $(SRNABENCH_LIBS)/sRNAbenchOutputDescription.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/sRNAbenchOutputDescription.txt 
	## END PIPELINE
	@echo -e "$(ts) SMRNAPIPELINE: END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: END\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "#END OF STATS from smallRNA-seq Pipeline for sample $(SAMPLE_ID) \n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
##
##