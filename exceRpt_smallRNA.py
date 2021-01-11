import os
import sys
import getpass
import shutil
import time
import pandas as pd

# RNA-seq pipeline usage:
help = '''

exceRpt smallRNA-seq pipeline

An example command to run the pipeline is:
    BASE=/Users/rk504/Pipelines/exceRpt_smallRNA.py
    PATH_OUT=/Users/rk504/tmp
    snakemake -s $BASE/snakemake_pipeline.py \
              --cores 2 \
              --directory=$PATH_OUT \
              --configfile="$BASE/TESTING/config_mac.yaml" 
              
Add the '-n -p' flags to do a dry run (and print the shell commands)

The two required input files are:
    - a yaml config file (passed to snakemake via the --configfile parameter)
    - a csv sample table (specified in the yaml config file via the "sample_table" : [path/tp/sampleTable.csv])

The output is written to a subdirectory of --directory=/path/to/output

An example of a valid config file is:
    "path_S3_in" : "s3://kitchen-mgh-data/Projects/TestData/RNAseq/Fastq"
    "path_S3_out" : "s3://kitchen-mgh-data/Projects/TestData/RNAseq/Processed"
    "sample_table" : "/Users/rk504/Pipelines/rnaseq-pipeline/TESTING/sampleTable.csv"
    "path_annotation" : "/Users/rk504/Dropbox_MGH/Annotations/Human"
    "data_origin": ["file","S3","sra"] (choose one)
    "path_tmp" : "/tmp/RNAseq"
    "rapid_run" : "True"
    "additional_args_salmon" : "--some other salmon args"

And the sample table must be formatted as such:
    sampleID,fastqName
    test,test_R1_001.fq.gz
    test,test_R2_001.fq.gz

In the sample table it is imperative that:
    - the fastqName is present at the S3 URI specified in 'path_S3_in' in the config file
'''


# hard-coded parameters
#exe_samtools = "/efs/bin/samtools-1.9/samtools"
#exe_bedtools = "/efs/bin/bedtools2/bin/bedtools"
#path_pipelineScripts = "/efs/Pipelines/RNAseq_pipeline/scripts"
job_name_sep = "-"


# Read parameters from the config file
path_in = config['path_in']
path_out = config['path_out']
#path_EFS_out = config['path_EFS_out']
path_ann = config["path_annotation"]
data_source = config["data_source"]
path_bin = config["path_bin"]
path_exceRpt = config["path_exceRpt"]

path_tmp = "/tmp"
if "path_tmp" in config:
    path_tmp = config["path_tmp"]

run_rapid = "False"
if "rapid_run" in config:
    run_rapid = config["rapid_run"]

additional_args_salmon = ""
if "additional_args_salmon" in config:
    additional_args_salmon = config["additional_args_salmon"]


exe_samtools = path_bin + "/samtools-1.9/samtools"
exe_bedtools = path_bin + "/bedtools2/bin/bedtools"
exe_star = path_bin + "/STAR/bin/STAR"

# Get a list of the samples and read numbers we need to process
samples = pd.read_csv(config["sample_table"], comment='#')
unique_sampleIDs = samples['sampleID'].unique()
unique_readNumbers = "R1"


# Get the system user
sys_username = getpass.getuser()


# Function to get the fastq filenames for all samples in the sample table
# -> returns a list of --include= statements for downloading from the S3 bucket via 'aws s3 sync'
def get_fastq(wildcards):
    tmp = samples.loc[(samples['sampleID'] == wildcards.sampleID)]
    fileNames = tmp['fastqName'].values
    return fileNames


# Primary rule - determine which sub-rules below are executed
if not run_rapid == "True":
    rule all:
        input:
            expand("{sampleID}/{sampleID}.{read}_fastqc.html",
                   sampleID=unique_sampleIDs, read=unique_readNumbers),
            expand("{sampleID}/trimmed/{sampleID}.{read}_fastqc.html",
                   sampleID=unique_sampleIDs, read=unique_readNumbers),
            #expand("{sampleID}/salmon/{sampleID}.coverage.csv.gz", sampleID=unique_sampleIDs),
            expand("{sampleID}/checkpoints/cleanup_and_sync.chk",
                   sampleID=unique_sampleIDs)
else:
    rule all:
        input:
            expand("{sampleID}/checkpoints/cleanup_and_sync.chk",
                   sampleID=unique_sampleIDs)


rule cleanup_and_sync:
    input:
        # "{sampleID}/checkpoints/sortBAM_pos.chk",
        # "{sampleID}/checkpoints/sortBAM_name.chk",
        "{sampleID}/checkpoints/calculate_stats.chk"
    output:
        "{sampleID}/checkpoints/cleanup_and_sync.chk"
    params:
        threads = 1,
        usethreads = 1,
        sampleID = "{sampleID}",
        runtime = "00:30:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"clean_up"
    # shell:
    #     '''
    #     rm -f {params.sampleID}/{params.sampleID}.R1.fastq.gz \
    #     && rm -f -r {params.sampleID}/rawFastq \
    #     && rm -f -r {params.sampleID}/trimmed \
    #     && rm -f {params.sampleID}/salmon/{params.sampleID}.bam \
    #     && aws s3 sync {params.sampleID} {path_S3_out}/{params.sampleID} --quiet \
    #     && rm -f {params.sampleID}/salmon/*.bam* \
    #     && touch {output} \
    #     && cp -r {params.sampleID} {path_EFS_out}/
    #     '''
    run:
        shell("$(COMPRESS_COMMAND)")
        shell("tar -cvz -C $(OUTPUT_DIR) -T $(OUTPUT_DIR)/$(SAMPLE_ID)_filesToCompress.txt -f $(OUTPUT_DIR)/$(SAMPLE_ID)_CORE_RESULTS_v$(EXCERPT_VERSION).tgz 2> /dev/null")


rule calculate_stats:
    input:
        a = "{sampleID}/checkpoints/process_alignments.chk",
        b = "{sampleID}/checkpoints/calculate_sequence_lengths_R1.chk",
        c = "{sampleID}/trimmed/{sampleID}.R1_fastqc.html"
    output:
        chk = "{sampleID}/checkpoints/calculate_stats.chk",
        stats = "{sampleID}/sampleStats.csv"
    params:
        threads = 1,
        usethreads = 1,
        sampleID = "{sampleID}",
        runtime = "00:10:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"calculate_stats"
    log:
        "{sampleID}/logs/calculate_stats.log"
    shell:
        '''
        echo "What,ReadCount" > {output.stats} \
        && cat {params.sampleID}/logs/trim_adapters.log | grep -w "^Input:" \
            | awk -F " " '{{print "Input,"$2}}' >> {output.stats} \
        && cat {params.sampleID}/logs/trim_adapters.log | grep -w "^Result:" \
            | awk -F " " '{{print "TrimmedFiltered,"$2}}' >> {output.stats} \
        && cat {params.sampleID}/salmon/aux_info/meta_info.json | grep -w "num_mapped" \
            | sed 's/,//g' | awk -F " " '{{print "TranscriptomeMapped,"$2}}' >> {output.stats} \
        && touch {output.chk}
        '''


rule process_alignments:
    input:
        a = "{sampleID}/checkpoints/map_RNA_genomeMapped.chk",
        b = "{sampleID}/checkpoints/map_RNA_genomeUnmapped.chk"
    output:
        "{sampleID}/checkpoints/process_alignments.chk"
    params:
        threads = 8,
        usethreads = 8,
        sampleID = "{sampleID}",
        runtime = "02:00:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"process_alignments"
    log:
        "{sampleID}/logs/process_alignments.log"
    shell:
        '''
        mkdir -p {path_tmp}/{params.sampleID} \
        && gunzip -c {path_ann}/annotation.gtf.gz > {path_tmp}/{params.sampleID}/annotation.gtf \
        && docker run -u `id -u {sys_username}` \
            -v $PWD/{params.sampleID}:/base -v {path_ann}:/ann -v {path_tmp}/{params.sampleID}:/scratch \
            rkitchen/salmon quant -i /ann/salmonIndex_31nt -g /scratch/annotation.gtf \
            --threads={params.usethreads} --seqBias --gcBias -l A --numBootstraps=20 \
            {additional_args_salmon} \
            --validateMappings --allowDovetail --writeMappings=/scratch/{params.sampleID}.sam \
            -o /base/salmon -r /base/trimmed/{params.sampleID}.R1.fastq.gz \
            >> {log} 2>&1 \
        && touch {params.sampleID}/checkpoints/process_alignments.chk \
        && {exe_samtools} view -bS --threads {params.threads} \
            {path_tmp}/{params.sampleID}/{params.sampleID}.sam \
            > {params.sampleID}/salmon/{params.sampleID}.bam \
        && rm -r {path_tmp}/{params.sampleID} \
        && touch {params.sampleID}/checkpoints/quantify_transcripts.chk \
        '''

rule map_RNA_genomeUnmapped:
    input:
        "{sampleID}/checkpoints/map_genome.chk"
    output:
        "{sampleID}/checkpoints/map_RNA_genomeUnmapped.chk"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        #read = "{read}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep + \
            "map_RNA_genomeUnmapped"
    log:
        "{sampleID}/logs/map_RNA_genomeUnmapped.log"
    shell:
        '''
        '''


rule map_RNA_genomeMapped:
    input:
        "{sampleID}/checkpoints/map_genome.chk"
    output:
        "{sampleID}/checkpoints/map_RNA_genomeMapped.chk"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        #read = "{read}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep + \
            "map_RNA_genomeMapped"
    log:
        "{sampleID}/logs/map_RNA_genomeMapped.log"
    shell:
        '''
        '''


rule map_genome:
    input:
        "{sampleID}/checkpoints/map_univec_and_rRNA.chk"
    output:
        "{sampleID}/checkpoints/map_genome.chk"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        #read = "{read}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep+"map_genome"
    log:
        "{sampleID}/logs/map_genome.log"
    shell:
        '''
        '''

rule map_rRNA_and_UniVec:
    input:
        "{sampleID}/checkpoints/trim_adapters_R1.chk"
    output:
        "{sampleID}/checkpoints/map_univec_and_rRNA.chk"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        #read = "{read}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep+"map_rRNA_and_UniVec"
    log:
        "{sampleID}/logs/map_rRNA_and_UniVec.log"
    shell:
        '''
        {exe_star} --runThreadN $(N_THREADS) \
          --outFileNamePrefix $(OUTPUT_DIR)/$(SAMPLE_ID)/filteringAlignments_UniVec_and_rRNA_ \
          --genomeDir $(DATABASE_PATH)/$(MAIN_ORGANISM_GENOME_ID)/STAR_INDEX_Univec_rRNA \
          --readFilesIn $(FILE_TO_INPUT_TO_UNIVEC_ALIGNMENT) --readFilesCommand "gunzip -c" \
          --outReadsUnmapped Fastx --parametersFiles $(DATABASE_PATH)/STAR_Parameters_Endogenous_smallRNA.in \
          $(STAR_ENDOGENOUS_DYNAMIC_PARAMS) >{log} 2>&1; \
        && {exe_samtools} view $(OUTPUT_DIR)/$(SAMPLE_ID)/filteringAlignments_UniVec_and_rRNA_Aligned.out.bam \
          | awk "{{print $3}}" | sort -k 2,2 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err \
          | uniq -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.trimmed.filtered.UniVec_and_rRNA.counts \
          2>>{log}; \
        && {exe_samtools} view $(OUTPUT_DIR)/$(SAMPLE_ID)/filteringAlignments_UniVec_and_rRNA_Aligned.out.bam \
          | awk "{{print $1}}" | sort 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err | uniq -c | wc -l \
          > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.trimmed.filtered.UniVec_and_rRNA.readCount \
          2>>{log}; \
        && gzip -c $(OUTPUT_DIR)/$(SAMPLE_ID)/filteringAlignments_UniVec_and_rRNA_Unmapped.out.mate1 \
          > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).clipped.trimmed.filtered.noUniVecOrRiboRNA.fastq.gz; \
        && rm $(OUTPUT_DIR)/$(SAMPLE_ID)/filteringAlignments_UniVec_and_rRNA_Unmapped.out.mate1
        '''

rule fastQC_trimmed:
    input:
        "{sampleID}/checkpoints/trim_adapters_{read}.chk"
    output:
        "{sampleID}/trimmed/{sampleID}.{read}_fastqc.html"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        read = "{read}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_{read}"+job_name_sep+"fastQC_trimmed"
    log:
        "{sampleID}/logs/fastQC_trimmed_{read}.log"
    shell:
        '''
        docker run -u `id -u {sys_username}` -v $PWD/{params.sampleID}:/base rkitchen/fastqc \
            -o /base/trimmed/ -t {params.usethreads} /base/trimmed/{params.sampleID}.{params.read}.fastq.gz \
            >> {log} 2>&1
        '''


rule calculate_sequence_lengths:
    input:
        "{sampleID}/checkpoints/trim_adapters_R1.chk"
    output:
        "{sampleID}/checkpoints/calculate_sequence_lengths_R1.chk"
    params:
        threads = 8,
        usethreads = 8,
        sampleID = "{sampleID}",
        runtime = "02:00:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"calculate_sequence_lengths"
    log:
        "{sampleID}/logs/calculate_sequence_lengths.log"
    shell:
        '''
        java -Xmx10G -jar {path_exceRpt}/exceRpt_Tools.jar GetSequenceLengths \
        $PWD/{params.sampleID}/trimmed/{params.sampleID}.fastq.gz \
        > $PWD/{params.sampleID}/trimmed/{params.sampleID}_insertSizes.txt
        '''

rule trim_adapters:
    input:
        R1 = "{sampleID}/checkpoints/combine_fastqs_R1.chk"
    output:
        R1 = "{sampleID}/checkpoints/trim_adapters_R1.chk"
    params:
        threads = 8,
        usethreads = 8,
        sampleID = "{sampleID}",
        runtime = "02:00:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"trim_adapters"
    log:
        "{sampleID}/logs/trim_adapters.log"
    shell:
        '''
        docker run --entrypoint bbduk.sh -u `id -u {sys_username}` -v $PWD/{params.sampleID}:/base rkitchen/bbtools \
            ref=/usr/local/bin/bbmap/resources/adapters.fa \
            in=/base/{params.sampleID}.R1.fastq.gz \
            out=/base/trimmed/{params.sampleID}.R1.fastq.gz \
            ktrim=r k=21 mink=11 tbo tpe hdist=2 minlen=31 threads={params.usethreads} \
            qtrim=r trimq=10 maq=10 \
            entropy=0.3 entropywindow=50 entropyk=5 \
            >> {log} 2>&1 \
        && touch {output.R1}
        '''


rule combine_fastqs:
    output:
        "{sampleID}/checkpoints/combine_fastqs_{read}.chk"
    params:
        sampleID = "{sampleID}",
        read = "{read}",
        files = get_fastq,
        threads = 1,
        usethreads = 1,
        runtime = "01:00:00",
        priority = 10,
        name = "{sampleID}_{read}"+job_name_sep+"combine_fastqs"
    # log:
    #    "{sampleID}/logs/combine_fastqs_{read}.log"
    run:
        shell("mkdir -p {params.sampleID}/logs")
        shell("mkdir -p {params.sampleID}/checkpoints")
        shell(
            "rm -f -r {params.sampleID}/{params.sampleID}.{params.read}.fastq.gz")
        for f in params.files:
            if data_origin == "S3":
                shell('aws s3 cp "{path_in}/' + f +
                      '" - >> {params.sampleID}/{params.sampleID}.{params.read}.fastq.gz')
            elif data_origin == "sra":
                shell('{exe_sra} --stdout "' + f +
                      '" >> {params.sampleID}/{params.sampleID}.{params.read}.fastq.gz')
            else:
                shell('cat "{path_in}/' + f +
                      '" >> {params.sampleID}/{params.sampleID}.{params.read}.fastq.gz')
        shell("touch {output}")
