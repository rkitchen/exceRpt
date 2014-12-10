#######################################################################################
##                                                                                   ##
## SmallRNA-seq pipeline - processes a single sequence file from a single sample     ##
##                                                                                   ##
## Author: Rob Kitchen (rob.kitchen@yale.edu)                                        ##
##                                                                                   ##
## Version 1.3.2 (2014-12-10)                                                        ##
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
## 2) Choose the main organism for smallRNA / genome alignment
##
## Human:
MAIN_ORGANISM           := hsa
#MAIN_ORGANISM_GENOME_ID := hg19
MAIN_ORGANISM_GENOME_ID := hg38
## Mouse:
#MAIN_ORGANISM           := mmu
#MAIN_ORGANISM_GENOME_ID := mm10


##
## 3) Select whether pipeline is run locally, should be 'true' unless this is the Genboree implementation!
##
LOCAL_EXECUTION := true


##
### If this is a local installation of the pipeline, be sure to also modify the parameters in steps 4, 5, and 6 below...
##
ifeq ($(LOCAL_EXECUTION),true)

	##
	## 4) Modify installation-specific variables
	##
	N_THREADS := 4
	MAX_RAM   := 12000000000
	## NB: The 'EXE_DIR' MUST be an ABSOLUTE PATH or sRNABench will fail!
	EXE_DIR   := /home/fas/gerstein/rrk24/scratch/bin/smallRNAPipeline
	
	
	##
	## 5) Choose optional analysis-specific options (or specify these at the command line)
	##
	ADAPTER_SEQ             := NULL
	TRNA_MAPPING            := on
	SNORNA_MAPPING          := on
	PIRNA_MAPPING           := on
	RFAM_MAPPING            := on
	GENCODE_MAPPING         := on
	MAP_PLANTS_VIRUSES      := on
	
	##
	## 6) Check that the paths to the required 3rd party executables work!
	##
	JAVA_EXE       := /usr/bin/java
	FASTX_EXE      := $(EXE_DIR)/fastx_0.0.13/fastx_clipper_33
	#BOWTIE1_PATH  := $(EXE_DIR)/bowtie-1.1.1
	#VIENNA_PATH   := $(EXE_DIR)/ViennaRNA_2.1.5/bin
	BOWTIE_EXE     := $(EXE_DIR)/bowtie2-2.1.0/bowtie2
	SAMTOOLS_EXE   := $(EXE_DIR)/samtools-0.1.18/samtools
	FASTQC_EXE     := $(JAVA_EXE) -classpath $(EXE_DIR)/FastQC_0.9.4
	SRATOOLS_EXE   := $(EXE_DIR)/sratoolkit.2.1.7-centos_linux64/fastq-dump
	THUNDER_EXE    := $(EXE_DIR)/Thunder.jar
	SRNABENCH_EXE  := $(EXE_DIR)/sRNAbench.jar
	SRNABENCH_LIBS := $(EXE_DIR)/sRNAbenchDB
	
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
	MAX_RAM := 12000000000

	FASTX_EXE := fastx_clipper_33
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
BOWTIE_OVERRIDE := \#2\#9\#v\#100

ifeq ($(MAIN_ORGANISM),hsa)  ## FOR HUMAN
	
	## Override the genome for adapter identification (saves having bt2 indexes for both hg19 and hg38)
	GENOME_ID_FOR_ADAPTER := hg38

	ifeq ($(MAIN_ORGANISM_GENOME_ID),hg19) ## hg19

		INDEX_GENCODE 		:= hg19_gencode18
		INDEX_TRNA			:= hg19_gencode18_tRNAs
		INDEX_PIRNA	 		:= hg19_piRNAs
		BOWTIE_INDEX_RRNA	:= $(SRNABENCH_LIBS)/customIndices/hg19_rRNA
		
	else ifeq ($(MAIN_ORGANISM_GENOME_ID),hg38) ## hg38
		
		INDEX_GENCODE 		:= hg38_gencode21
		INDEX_TRNA			:= hg38_gencode21_tRNAs
		INDEX_PIRNA	 		:= hg19_piRNAs
		BOWTIE_INDEX_RRNA	:= $(SRNABENCH_LIBS)/customIndices/hg38_rRNA
		
	endif

	GENCODE_LIBS := libs=$(INDEX_GENCODE)$(BOWTIE_OVERRIDE)
	TRNA_LIBS	 := libs=$(INDEX_TRNA)$(BOWTIE_OVERRIDE)
	PIRNA_LIBS   := libs=$(INDEX_PIRNA)$(BOWTIE_OVERRIDE)

else ifeq ($(MAIN_ORGANISM),mmu)  ## FOR MOUSE
	
	GENOME_ID_FOR_ADAPTER := mm10

	INDEX_GENCODE 		:= mm10_gencodeM4
	INDEX_TRNA			:= mm10_gencodeM4_tRNAs
	INDEX_PIRNA	 		:= mm10_piRNAs
	BOWTIE_INDEX_RRNA	:= $(SRNABENCH_LIBS)/customIndices/mm10_rRNA

	GENCODE_LIBS := libs=$(INDEX_GENCODE)$(BOWTIE_OVERRIDE)
	TRNA_LIBS	 := libs=$(INDEX_TRNA)$(BOWTIE_OVERRIDE)
	PIRNA_LIBS   := libs=$(INDEX_PIRNA)$(BOWTIE_OVERRIDE)
	
endif


##
## Turn off a smallRNA library if they are not selected by the user
##
ifneq ($(TRNA_MAPPING),on)
	TRNA_LIBS :=
endif
ifneq ($(SNORNA_MAPPING),on)
	SNORNA_LIBS :=
endif
ifneq ($(PIRNA_MAPPING),on)
	PIRNA_LIBS :=
endif
ifneq ($(RFAM_MAPPING),on)
	RFAM_LIBS :=
endif
ifneq ($(GENCODE_MAPPING),on)
	GENCODE_LIBS :=
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
OTHER_LIBRARIES := $(TRNA_LIBS) $(SNORNA_LIBS) $(PIRNA_LIBS) $(RFAM_LIBS) $(GENCODE_LIBS)

## Map reads to plant and virus miRNAs
ifeq ($(MAP_PLANTS_VIRUSES),on)
	#PROCESS_SAMPLE_REQFILE := PlantAndVirus/mature_sense_nonRed.grouped
	PROCESS_SAMPLE_REQFILE := PlantAndVirus/mature_sense.grouped
else
	PROCESS_SAMPLE_REQFILE := readsNotAssigned.fa
endif


##
## List of plant and virus species IDs to which to map reads that do not map to the genome of the primary organism
##
PLANT_LIST := cre:cln:pab:pde:pta:ppt:smo:pgi:cca:han:har:hci:hex:hpa:hpe:htu:aly:ath:bna:bol:bra:cpa:cme:hbr:mes:rco:aau:ahy:amg:gma:gso:lja:mtr:pvu:vun:ama:dpr:rgl:ssl:lus:gar:ghb:ghr:gra:tcc:aqc:bcy:bgy:mdm:ppe:ccl:crt:csi:ctr:peu:ptc:nta:sly:stu:vvi:ata:bdi:egu:far:hvu:osa:sbi:sof:ssp:tae:ttu:zma
VIRUS_LIST := bhv1:bkv:blv:bpcv1:bpcv2:dev:ebv:hbv:hcmv:hhv6b:hiv1:hsv1:hsv2:hvsa:hvt:iltv:jcv:kshv:mcmv:mcv:mdv1:mdv2:mghv:prv:rlcv:rrv:sv40

## Parameters to use for the bowtie mapping of calibrator oligos and rRNAs
#BOWTIE2_MAPPING_PARAMS := -D 15 -R 2 -N 1 -L 16 -i S,1,1.15
BOWTIE2_MAPPING_PARAMS_CALIBRATOR := -D 15 -R 2 -N 1 -L 16 -i S,1,0
BOWTIE2_MAPPING_PARAMS_RRNA := -D 15 -R 2 -N 1 -L 19 -i S,1,0

#################################################

USEAGE := 
ifeq ($(INPUT_FILE_ID),NULL)
  USEAGE := "make -f smallRNA_pipeline INPUT_FILE_PATH=[required: absolute/path/to/input/.fa|.fq|.sra] N_THREADS=[required: number of threads] OUTPUT_DIR=<required: absolute/path/to/output> INPUT_FILE_ID=[required: samplename] ADAPTER_SEQ=[optional: will guess sequence if not provided here; none, if already clipped input] MAIN_ORGANISM=[optional: defaults to 'hsa'] MAIN_ORGANISM_GENOME_ID=[optional: defaults to 'hg19'] CALIBRATOR_LIBRARY=[optional: path/to/bowtie/index/containing/calibrator/sequences] TRNA_MAPPING=[optional: TRUE|FALSE, default is TRUE] SNORNA_MAPPING=[optional: TRUE|FALSE, default is TRUE] PIRNA_MAPPING=[optional: TRUE|FALSE, default is TRUE] MAP_RFAM=[optional: TRUE|FALSE, default is TRUE] MAP_PLANTS_VIRUSES=[optional: TRUE|FALSE, default is TRUE]"
endif



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
COMMAND_CLIP_ADAPTER ?= $(COMMAND_CONVERT_SRA) | $(FASTX_EXE) -a $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq) -l 10 -v -M 7 -z -o $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err




##
## Logic block for removing rRNAs and [optionally] calibrator sequences that may have been spiked into the sample
##
ifeq ($(CALIBRATOR_LIBRARY),NULL)

	LOGENTRY_MAP_CALIBRATOR_1 := No calibrator sequences\n 
	LOGENTRY_MAP_CALIBRATOR_2 := Moving on to UniVec and rRNA sequences\n
	COMMAND_COUNT_CALIBRATOR := echo -e "calibrator\tNA" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	COMMAND_MAP_CALIBRATOR := 
	
	FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT := $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz
	#COMMAND_MAP_RRNAS := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noRiboRNA.fastq.gz -x $(BOWTIE_INDEX_RRNA) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' > tmp.sam; \
	#cat tmp.sam | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
	#cat tmp.sam | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.bam; \
	#rm tmp.sam
else

	COMMAND_COUNT_CALIBRATOR := cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.calibratormapped.counts | awk '{sum+=$$1} END {print "calibrator\t"sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	COMMAND_MAP_CALIBRATOR := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_CALIBRATOR) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noCalibrator.fastq.gz -x $(CALIBRATOR_LIBRARY) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | tee $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.calibratormapped.bam | $(SAMTOOLS_EXE) view - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.calibratormapped.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	LOGENTRY_MAP_CALIBRATOR_1 := $(ts) SMRNAPIPELINE: Mapping reads to calibrator sequences using bowtie:\n
	LOGENTRY_MAP_CALIBRATOR_2 := $(ts) SMRNAPIPELINE: Finished mapping to the calibrators\n
	
	FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT := $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noCalibrator.fastq.gz
	#COMMAND_MAP_RRNAS := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noRiboRNA.fastq.gz -x $(BOWTIE_INDEX_RRNA) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noCalibrator.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | tee $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.bam | $(SAMTOOLS_EXE) view - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
endif


##
## Bowtie2 command to align reads to the UniVec contaminant sequence database
##
COMMAND_MAP_UNIVEC := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noUniVecContaminants.fastq.gz -x $(BOWTIE_INDEX_UNIVEC) -U $(FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' > tmp.sam; \
cat tmp.sam | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.uniVecContaminants.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat tmp.sam | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.uniVecContaminantMapped.bam; \
rm tmp.sam


##
## Bowtie2 command to align reads to the rRNA sequences
##
COMMAND_MAP_RRNAS := $(BOWTIE_EXE) -p $(N_THREADS) $(BOWTIE2_MAPPING_PARAMS_RRNA) --un-gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noRiboRNA.fastq.gz -x $(BOWTIE_INDEX_RRNA) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noUniVecContaminants.fastq.gz 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '$$2 != 4 {print $$0}' > tmp.sam; \
cat tmp.sam | awk '{print $$3}' | sort -k 2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq --count > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.counts 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err; \
cat tmp.sam | $(SAMTOOLS_EXE) view -Sb - 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).log > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.bam; \
rm tmp.sam



## Parameters for the endogenous-exRNA mapping
FIXED_PARAMS_MAIN := -Xmx10G -jar $(SRNABENCH_EXE) dbPath=$(SRNABENCH_LIBS) p=$(N_THREADS) chunkmbs=2000 microRNA=$(MAIN_ORGANISM) species=$(MAIN_ORGANISM_GENOME_ID) plotMiR=true plotLibs=true predict=false $(OTHER_LIBRARIES) writeGenomeDist=true noMM=1 maxReadLength=75

## Parameters for the exogenous-exRNA mapping
FIXED_PARAMS_PLANT_VIRUS := -Xmx10G -jar $(SRNABENCH_EXE) dbPath=$(SRNABENCH_LIBS) p=$(N_THREADS) chunkmbs=2000 microRNA=$(PLANT_LIST):$(VIRUS_LIST) plotMiR=true predict=false noMM=1










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
####  vvv Sub-targets to do the read-preprocessing, calibrator mapping, rRNA mapping, en-exRNA mapping, and ex-exRNA mapping vvv
###
##


## BEGIN PIPELINE
##
## Make results directory
##
$(OUTPUT_DIR)/$(SAMPLE_ID): 
	#$(EXPORT_CMD)
	@echo -e "$(USEAGE)"
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "$(ts) SMRNAPIPELINE: BEGIN smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" > $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: BEGIN \n" > $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "#STATS from smallRNA-seq Pipeline for sample $(SAMPLE_ID)" > $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	@echo -e "$(ts) SMRNAPIPELINE: Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log



##
## Write adapter sequence
##
#$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq: $(OUTPUT_DIR)/$(SAMPLE_ID)
$(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat: $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Processing adapter sequence:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_WRITE_ADAPTER_SEQ)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_WRITE_ADAPTER_SEQ)
	@echo -e "$(ts) SMRNAPIPELINE: Progress_1_FoundAdapter" > $(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat


##
## CLIP 3' adapter sequence
##
#$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/Progress_1_FoundAdapter.dat
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Adapter sequence: $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).adapterSeq)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(LOGENTRY_WRITE_ADAPTER)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_CLIP_ADAPTER)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_CLIP_ADAPTER)
	@echo -e "$(ts) SMRNAPIPELINE: Finished removing adapters\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Assess Read-lengths after clipping
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.readLengths.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz
	@echo -e "======================" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Calculating length distribution of clipped reads:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	gunzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq
	@echo -e "$(ts) SMRNAPIPELINE: $(JAVA_EXE) -Xmx1G -jar $(THUNDER_EXE) GetSequenceLengths $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.readLengths.txt 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(JAVA_EXE) -Xmx1G -jar $(THUNDER_EXE) GetSequenceLengths $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.readLengths.txt 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq
	@echo -e "$(ts) SMRNAPIPELINE: Finished calculating read-lengths\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log



##
## MAP to external bowtie (calibrator?) library and to the rRNA sequences
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noRiboRNA.fastq.gz: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.fastq.gz $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.readLengths.txt
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(LOGENTRY_MAP_CALIBRATOR_1)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_MAP_CALIBRATOR)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_CALIBRATOR)
	@echo -e "$(ts) SMRNAPIPELINE: $(LOGENTRY_MAP_CALIBRATOR_2)" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to contaminant sequences in UniVec using Bowtie2:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_MAP_UNIVEC)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_UNIVEC)
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to the UniVec contaminant DB\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to ribosomal RNA sequences using Bowtie2:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(COMMAND_MAP_RRNAS)\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	$(COMMAND_MAP_RRNAS) 
	$(SAMTOOLS_EXE) sort $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.bam $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.sorted
	$(SAMTOOLS_EXE) index $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.sorted.bam
	rm $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.bam
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to the rRNAs\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Map reads to the main genome of interest
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/readsNotAssigned.fa: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noRiboRNA.fastq.gz
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to smallRNAs of primary organism:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(JAVA_EXE) $(FIXED_PARAMS_MAIN) input=$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noRiboRNA.fastq.gz output=$(OUTPUT_DIR)/$(SAMPLE_ID) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 
	$(JAVA_EXE) $(FIXED_PARAMS_MAIN) input=$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.noRiboRNA.fastq.gz output=$(OUTPUT_DIR)/$(SAMPLE_ID) >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to the small-RNAs of the primary organism\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: grep "chrUn_gl000220" $(OUTPUT_DIR)/$(SAMPLE_ID)/genome.txt | awk '{print $$1}' | sed 's/#/ /g' | awk '{ sum += $$2 } END { print "Number of reads mapped to chrUn_gl000220 = "sum }' >> $(OUTPUT_DIR)/$(SAMPLE_ID).log\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	grep "chrUn_gl000220" $(OUTPUT_DIR)/$(SAMPLE_ID)/genome.txt | awk '{print $$1}' | sed 's/#/ /g' | awk '{ sum += $$2 } END { print "Number of reads mapped to chrUn_gl000220 = "sum }' >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## Use the unmapped reads and search against all plant and viral miRNAs in miRBase
##
#$(OUTPUT_DIR)/$(SAMPLE_ID)/PlantAndVirus/mature_sense_nonRed.grouped: $(OUTPUT_DIR)/$(SAMPLE_ID)/readsNotAssigned.fa
$(OUTPUT_DIR)/$(SAMPLE_ID)/PlantAndVirus/mature_sense.grouped: $(OUTPUT_DIR)/$(SAMPLE_ID)/readsNotAssigned.fa
	@echo -e "======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: Mapping reads to smallRNAs of plants and viruses:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: $(JAVA_EXE) $(FIXED_PARAMS_PLANT_VIRUS) input=$(OUTPUT_DIR)/$(SAMPLE_ID)/readsNotDetected.fa output=$(OUTPUT_DIR)/$(SAMPLE_ID)/PlantAndVirus >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	#$(JAVA_EXE) $(FIXED_PARAMS_PLANT_VIRUS) input=$(OUTPUT_DIR)/$(SAMPLE_ID)/readsNotDetected.fa output=$(OUTPUT_DIR)/$(SAMPLE_ID)/PlantAndVirus >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(JAVA_EXE) $(FIXED_PARAMS_PLANT_VIRUS) input=$(OUTPUT_DIR)/$(SAMPLE_ID)/readsNotAssigned.fa output=$(OUTPUT_DIR)/$(SAMPLE_ID)/PlantAndVirus >> $(OUTPUT_DIR)/$(SAMPLE_ID).log 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "$(ts) SMRNAPIPELINE: Finished mapping to plant and virus small-RNAs\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

##
## Summary stats
##
processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(PROCESS_SAMPLE_REQFILE)
#processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/PlantAndVirus/mature_sense_nonRed.grouped
	@echo -e "Stage\tReadCount" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

	## Count reads input to adapter clipping
	grep "Input: " $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "input\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	
	## Count reads output following adapter clipping
	grep "Output: " $(OUTPUT_DIR)/$(SAMPLE_ID).log | awk '{print "clipped\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

	## Count calibrator oligo reads
	$(COMMAND_COUNT_CALIBRATOR)
    
	## Count UniVec contaminant reads
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.uniVecContaminants.counts | awk '{sum+=$$1} END {print "UniVec_contaminants\t"sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

	## Count rRNA reads
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.rRNAmapped.counts | awk '{sum+=$$1} END {print "rRNA\t"sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

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

	## Assigned non-redundantly to annotated snoRNAs
	#cat $(OUTPUT_DIR)/$(SAMPLE_ID)/snorna_sense.grouped | awk '{sum+=$$4} END {printf "snoRNA_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	#cat $(OUTPUT_DIR)/$(SAMPLE_ID)/snorna_antisense.grouped | awk '{sum+=$$4} END {printf "snoRNA_antisense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

	## Assigned non-redundantly to annotated transcripts in Gencode
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_GENCODE)_sense.grouped | awk '{sum+=$$4} END {printf "Gencode_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(INDEX_GENCODE)_antisense.grouped | awk '{sum+=$$4} END {printf "Gencode_antisense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

	## Assigned non-redundantly to annotated plant/viral miRNAs
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/PlantAndVirus/mature_sense.grouped | awk '{sum+=$$4} END {printf "miRNA_plantVirus_sense\t%.0f\n",sum}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

	## Copy Output descriptions file
	cp $(SRNABENCH_LIBS)/sRNAbenchOutputDescription.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/sRNAbenchOutputDescription.txt 
	
## END PIPELINE
	@echo -e "$(ts) SMRNAPIPELINE: END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log
	@echo -e "$(ts) SMRNAPIPELINE: END\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	@echo -e "#END OF STATS from smallRNA-seq Pipeline for sample $(SAMPLE_ID) \n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats

##
##
