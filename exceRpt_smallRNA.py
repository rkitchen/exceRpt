import os
import sys
import getpass
import shutil
import time
import pandas as pd

version_exceRpt = "5.1"

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
# exe_samtools = "/efs/bin/samtools-1.9/samtools"
# exe_bedtools = "/efs/bin/bedtools2/bin/bedtools"
# path_pipelineScripts = "/efs/Pipelines/RNAseq_pipeline/scripts"
job_name_sep = "-"


# Read parameters from the config file
path_in = config['path_in']
path_out = config['path_out']
# path_EFS_out = config['path_EFS_out']
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

endogenous_lib_priority = "miRNA,tRNA,piRNA,gencode,circRNA"


minInsertSize = config["minInsertSize"]
min_adapter_bases_3p = 7
trim_bases_5p = 0
trim_bases_3p = 0


exe_samtools = path_bin + "/samtools-1.9/samtools"
exe_bedtools = path_bin + "/bedtools2/bin/bedtools"
# exe_star = path_bin + "/STAR/STAR-2.7.7a/bin/Linux_x86_64/STAR"
exe_star = path_bin + "/STAR/STAR-2.7.1a/bin/Linux_x86_64/STAR"


# Get a list of the samples and read numbers we need to process
samples = pd.read_csv(config["sample_table"], comment='#')
unique_sampleIDs = samples['sampleID'].unique()
unique_readNumbers = "R1"


params_star_endogenous = ' --alignEndsType ' + config["star_alignEndsType"] + ' \
        --outFilterMatchNmin ' + config["minInsertSize"] + ' \
        --outFilterMatchNminOverLread ' + config["star_outFilterMatchNminOverLread"] + ' \
        --outFilterMismatchNmax ' + config["star_outFilterMismatchNmax"] + ' \
        --outFilterMismatchNoverLmax ' + config["star_outFilterMismatchNoverLmax"] + ' '

params_star_exogenous = '--outSAMtype BAM Unsorted --outSAMattributes Standard --alignEndsType EndToEnd \
        --outFilterMatchNmin ' + config["minInsertSize"] + ' \
        --outFilterMatchNminOverLread ' + config["star_outFilterMatchNminOverLread"] + ' \
        --outFilterMismatchNmax ' + config["star_outFilterMismatchNmax"] + ' \
        --outFilterMismatchNoverLmax ' + config["star_outFilterMismatchNoverLmax"] + ' '


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
            # expand("{sampleID}/salmon/{sampleID}.coverage.csv.gz", sampleID=unique_sampleIDs),
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
        shell('''ls -lh {params.sampleID} | awk '{{print $9}}' \
                | grep "readCounts_\|.insertSizes.txt\|_fastqc.zip\|.counts\|.CIGARstats.txt\|.coverage.txt" \
                | awk '{{print "{params.sampleID}/"$1}}' \
                > {params.sampleID}/{params.sampleID}_filesToCompress.txt;''')
        shell("tar -cvz -C . -T {params.sampleID}/{params.sampleID}_filesToCompress.txt \
                -f {params.sampleID}_CORE_RESULTS_v{version_exceRpt}.tgz 2> /dev/null")
        shell('touch {output}')


rule calculate_stats:
    input:
        a = "{params.sampleID}/endogenousAlignments_Accepted.txt.gz",
        b = "{sampleID}/checkpoints/calculate_sequence_lengths.chk",
        c = "{sampleID}/{sampleID}_trimmed_filtered_fastqc.html"
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
    run:
        shell('''
              {exe_samtools} view {params.sampleID}/filteringAlignments_UniVec_and_rRNA_Aligned.out.bam \
              | awk '{{print $3}}' | sort -k 2,2 2>>{log} \
              | uniq -c | sort -nrk 1,1 \
              > {params.sampleID}/{params.sampleID}.clipped.trimmed.filtered.UniVec_and_rRNA.counts \
              2>{log}
        ''')
        shell('''
              {exe_samtools} view {params.sampleID}/filteringAlignments_UniVec_and_rRNA_Aligned.out.bam \
              | awk '{{print $1}}'| sort -nrk 1,1 2>>{log} | uniq -c | wc -l \
              > {params.sampleID}/{params.sampleID}.clipped.trimmed.filtered.UniVec_and_rRNA.readCount \
              2>>{log}
        ''')
        shell('''
              cat {params.sampleID}/readCounts_miRNAmature_sense.txt | awk '{{SUM+=$4}}END{{printf "miRNA_sense\t%.0f\n",SUM}}' >> {output.stats} \
              cat {params.sampleID}/readCounts_miRNAmature_antisense.txt | awk '{{SUM+=$4}}END{{printf "miRNA_antisense\t%.0f\n",SUM}}' >> {output.stats}
        ''')
        # Count reads not mapping to the genome or to the libraries
        shell('''
              gunzip -c {params.sampleID}/endogenousAlignments_genomeUnmapped_transcriptome_Unmapped.fastq.gz \
              | wc - l | awk '{{print "not_mapped_to_genome_or_libs\t"($1/4)}}' >> {output.stats}
              ''')
        shell('touch {output.chk}')

        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_miRNAprecursor_sense.txt | awk '{SUM+=$$4}END{printf "miRNAprecursor_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_miRNAprecursor_antisense.txt | awk '{SUM+=$$4}END{printf "miRNAprecursor_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_tRNA_sense.txt | awk '{SUM+=$$4}END{printf "tRNA_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_tRNA_antisense.txt | awk '{SUM+=$$4}END{printf "tRNA_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_piRNA_sense.txt | awk '{SUM+=$$4}END{printf "piRNA_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_piRNA_antisense.txt | awk '{SUM+=$$4}END{printf "piRNA_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_gencode_sense.txt | awk '{SUM+=$$4}END{printf "gencode_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_gencode_antisense.txt | awk '{SUM+=$$4}END{printf "gencode_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_circRNA_sense.txt | awk '{SUM+=$$4}END{printf "circularRNA_sense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats
        # cat $(OUTPUT_DIR) /$(SAMPLE_ID)/readCounts_circRNA_antisense.txt | awk '{SUM+=$$4}END{printf "circularRNA_antisense\t%.0f\n",SUM}' >> $(OUTPUT_DIR) /$(SAMPLE_ID).stats


rule process_alignments:
    input:
        a = "{sampleID}/checkpoints/map_RNA_genomeMapped.chk",
        b = "{sampleID}/checkpoints/map_RNA_genomeUnmapped.chk"
    output:
        chk = "{sampleID}/checkpoints/process_alignments.chk",
        alignments = "{sampleID}/endogenousAlignments_Accepted.txt.gz"
    params:
        threads = 1,
        usethreads = 1,
        javaMem = "10G",
        sampleID = "{sampleID}",
        runtime = "02:00:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"process_alignments"
    log:
        "{sampleID}/logs/process_alignments.log"
    run:
        shell('''
             java -Xmx{params.javaMem} -jar {path_exceRpt}/exceRpt_Tools.jar ProcessEndogenousAlignments \
             --libPriority {endogenous_lib_priority} \
             --genomeMappedReads {params.sampleID}/endogenousAlignments_genomeMapped_transcriptome_Aligned.out.bam \
             --transcriptomeMappedReads {params.sampleID}/endogenousAlignments_genomeUnmapped_transcriptome_Aligned.out.bam \
             --hairpin2genome {path_ann}/miRNA_precursor2genome.sam \
             --mature2hairpin {path_ann}/miRNA_mature2precursor.sam \
             --dict {params.sampleID}/endogenousAlignments_Accepted.dict \
             2>{log} | sort -k 2,2 -k 1,1 \
             > {params.sampleID}/endogenousAlignments_Accepted.txt
        ''')
        shell('''
              java -Xmx{params.javaMem} -jar {path_exceRpt}/exceRpt_Tools.jar QuantifyEndogenousAlignments \
              --dict {params.sampleID}/endogenousAlignments_Accepted.dict \
              --acceptedAlignments {params.sampleID}/endogenousAlignments_Accepted.txt \
              --outputPath {params.sampleID} 2>{log}
            ''')
        shell('gzip {params.sampleID}/endogenousAlignments_Accepted.txt')
        shell('touch {output}')


rule map_RNA_genomeUnmapped:
    input:
        "{sampleID}/endogenousAlignments_genome_Unmapped.fastq.gz"
    output:
        reads = "{sampleID}/endogenousAlignments_genomeUnmapped_transcriptome_Unmapped.fastq.gz",
        chk = "{sampleID}/checkpoints/map_RNA_genomeUnmapped.chk"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep + \
            "map_RNA_genomeUnmapped"
    log:
        "{sampleID}/logs/map_RNA_genomeUnmapped.log"
    run:
        shell('''
        {exe_star} --runThreadN {params.usethreads} \
        --outFileNamePrefix {params.sampleID}/endogenousAlignments_genomeUnmapped_transcriptome_ \
        --readFilesIn {input} \
        --outReadsUnmapped Fastx --genomeDir {path_ann}/STAR_INDEX_transcriptome \
        --parametersFiles {path_ann}/../STAR_Parameters_Endogenous_smallRNA.in \
        {params_star_endogenous} --readFilesCommand "gunzip -c" >{log} 2>&1
        ''')
        shell('mv {params.sampleID}/endogenousAlignments_genomeUnmapped_transcriptome_Unmapped.out.mate1 \
              {params.sampleID}/endogenousAlignments_genomeUnmapped_transcriptome_Unmapped.fastq')
        shell(
            'gzip {params.sampleID}/endogenousAlignments_genomeUnmapped_transcriptome_Unmapped.fastq')
        shell('touch {output.chk}')


rule map_RNA_genomeMapped:
    input:
        "{sampleID}/endogenousAlignments_genome_Aligned.out.bam"
    output:
        reads = "{sampleID}/endogenousAlignments_genomeMapped_transcriptome_Unmapped.fastq.gz",
        chk = "{sampleID}/checkpoints/map_RNA_genomeMapped.chk"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep +
        "map_RNA_genomeMapped"
    log:
        "{sampleID}/logs/map_RNA_genomeMapped.log"
    run:
        shell('''
            {exe_samtools} fastq {params.sampleID}/endogenousAlignments_genome_Aligned.out.bam \
            > {params.sampleID}/endogenousAlignments_genome_Mapped.out.mate1 2>{log}
        ''')
        shell('''
        {exe_star} --runThreadN {params.usethreads} \
        --outFileNamePrefix {params.sampleID}/endogenousAlignments_genomeMapped_transcriptome_ \
        --readFilesIn {params.sampleID}/endogenousAlignments_genome_Mapped.out.mate1 \
        --genomeDir {path_ann}/STAR_INDEX_transcriptome \
        --parametersFiles {path_ann}/../STAR_Parameters_Endogenous_smallRNA.in \
        {params_star_endogenous} --readFilesCommand - >>{log} 2>>{log}
        ''')
        shell('rm {params.sampleID}/endogenousAlignments_genome_Mapped.out.mate1')
        shell('mv {params.sampleID}/endogenousAlignments_genomeMapped_transcriptome_Unmapped.out.mate1 \
              {params.sampleID}/endogenousAlignments_genomeMapped_transcriptome_Unmapped.fastq')
        shell(
            'gzip {params.sampleID}/endogenousAlignments_genomeMapped_transcriptome_Unmapped.fastq')
        shell('touch {output.chk}')


rule map_genome:
    input:
        "{sampleID}/checkpoints/map_univec_and_rRNA.chk"
    output:
        bam = "{sampleID}/endogenousAlignments_genome_Aligned.out.bam",
        fastq_unmapped = "{sampleID}/endogenousAlignments_genome_Unmapped.fastq.gz"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep+"map_genome"
    log:
        "{sampleID}/logs/map_genome.log"
    run:
        shell('''
            {exe_star} --runThreadN {params.usethreads} \
            --outFileNamePrefix {params.sampleID}/endogenousAlignments_genome_ \
            --genomeDir {path_ann}/STAR_INDEX_genome \
            --readFilesIn {params.sampleID}/{params.sampleID}.clipped.trimmed.filtered.noUniVecOrRiboRNA.fastq.gz \
            --readFilesCommand "gunzip -c" --outReadsUnmapped Fastx \
            --parametersFiles {path_ann}/../STAR_Parameters_Endogenous_smallRNA.in \
            {params_star_endogenous} >{log} 2>&1
        ''')
        shell('mv {params.sampleID}/endogenousAlignments_genome_Unmapped.out.mate1 \
              {params.sampleID}/endogenousAlignments_genome_Unmapped.fastq')
        shell(
            'gzip {params.sampleID}/endogenousAlignments_genome_Unmapped.fastq')
        shell('touch {output}')


rule map_rRNA_and_UniVec:
    input:
        "{sampleID}/{sampleID}_trimmed_filtered.fastq.gz"
    output:
        "{sampleID}/checkpoints/map_univec_and_rRNA.chk"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}_"+job_name_sep+"map_rRNA_and_UniVec"
    log:
        "{sampleID}/logs/map_rRNA_and_UniVec.log"
    run:
        shell('''
        {exe_star} --runThreadN {params.usethreads} \
          --outFileNamePrefix {params.sampleID}/filteringAlignments_UniVec_and_rRNA_ \
          --genomeDir {path_ann}/STAR_INDEX_Univec_rRNA \
          --readFilesIn {input} --readFilesCommand "gunzip -c" \
          --outReadsUnmapped Fastx --parametersFiles {path_ann}/../STAR_Parameters_Endogenous_smallRNA.in \
          {params_star_endogenous} >{log} 2>&1
        ''')
        shell('''gzip -c {params.sampleID}/filteringAlignments_UniVec_and_rRNA_Unmapped.out.mate1 \
          > {params.sampleID}/{params.sampleID}.clipped.trimmed.filtered.noUniVecOrRiboRNA.fastq.gz''')
        shell(
            'rm {params.sampleID}/filteringAlignments_UniVec_and_rRNA_Unmapped.out.mate1')
        shell('touch {output}')


rule fastQC_trimmed:
    input:
        # "{sampleID}/checkpoints/trim_adapters.chk"
        "{sampleID}/{sampleID}_trimmed_filtered.fastq.gz"
    output:
        "{sampleID}/{sampleID}_trimmed_filtered_fastqc.html"
    params:
        threads = 4,
        usethreads = 4,
        sampleID = "{sampleID}",
        runtime = "01:00:00",
        priority = 1,
        name = "{sampleID}"+job_name_sep+"fastQC_trimmed"
    log:
        "{sampleID}/logs/fastQC_trimmed.log"
    shell:
        '''
        docker run -u `id -u {sys_username}` -v $PWD/{params.sampleID}:/base rkitchen/fastqc \
            -o /base -t {params.usethreads} /base/{params.sampleID}_trimmed_filtered.fastq.gz \
            >> {log} 2>&1
        '''


rule calculate_sequence_lengths:
    input:
        # "{sampleID}/checkpoints/trim_adapters.chk"
        "{sampleID}/{sampleID}_trimmed_filtered.fastq.gz"
    output:
        "{sampleID}/checkpoints/calculate_sequence_lengths.chk"
    params:
        threads = 8,
        usethreads = 8,
        sampleID = "{sampleID}",
        runtime = "02:00:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"calculate_sequence_lengths"
    log:
        "{sampleID}/logs/calculate_sequence_lengths.log"
    run:
        shell('''
        java -Xmx10G -jar {path_exceRpt}/exceRpt_Tools.jar GetSequenceLengths \
          {input} > {input}_insertSizes.txt
        ''')
        shell('touch {output}')


rule trim_adapters:
    input:
        R1 = "{sampleID}/checkpoints/combine_fastqs.chk"
    output:
        fastq = "{sampleID}/{sampleID}_trimmed_filtered.fastq.gz",
        chk = "{sampleID}/checkpoints/trim_adapters.chk"
    params:
        threads = 8,
        usethreads = 8,
        sampleID = "{sampleID}",
        runtime = "02:00:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"trim_adapters"
    log:
        "{sampleID}/logs/trim_adapters.log"
    run:
        shell('''
              docker run --entrypoint bbduk.sh -u `id -u {sys_username}` \
              -v $PWD/{params.sampleID}:/base rkitchen/bbtools \
              ref=/usr/local/bin/bbmap/resources/adapters.fa \
              in=/base/{params.sampleID}.fastq.gz \
              out=/base/{params.sampleID}_trimmed.fastq.gz \
              ktrim=r k=10 mink={min_adapter_bases_3p} trimpolya=7 hdist=0 minlen={minInsertSize} threads={params.usethreads} \
              forcetrimleft={trim_bases_5p} forcetrimright2={trim_bases_3p} \
              stats=/base/{params.sampleID}_filterStats_adapter.txt \
              >{log} 2>{log}
              ''')
        shell('''
              docker run --entrypoint bbduk.sh -u `id -u {sys_username}` \
              -v $PWD/{params.sampleID}:/base -v {path_ann}/..:/ann rkitchen/bbtools \
              ref=/ann/phiX.fa \
              in=/base/{params.sampleID}_trimmed.fastq.gz \
              out=/base/{params.sampleID}_trimmed_filtered.fastq.gz \
              qtrim=r trimq=10 maq=10 \
              entropy=0.3 entropywindow=50 entropyk=5 \
              k=20 hdist=1 minlen={minInsertSize} threads={params.usethreads} \
              stats=/base/{params.sampleID}_filterStats_phiXandQuality.txt \
              >>{log} 2>>{log}
              ''')
        shell('touch {output.chk}')


rule combine_fastqs:
    output:
        "{sampleID}/checkpoints/combine_fastqs.chk"
    params:
        sampleID = "{sampleID}",
        files = get_fastq,
        threads = 1,
        usethreads = 1,
        runtime = "01:00:00",
        priority = 10,
        name = "{sampleID}"+job_name_sep+"combine_fastqs"
    run:
        shell("mkdir -p {params.sampleID}/logs")
        shell("mkdir -p {params.sampleID}/checkpoints")
        shell(
            "rm -f -r {params.sampleID}/{params.sampleID}.fastq.gz")
        for f in params.files:
            if data_source == "S3":
                shell('aws s3 cp "{path_in}/' + f +
                      '" - >> {params.sampleID}/{params.sampleID}.fastq.gz')
            elif data_source == "sra":
                shell('{exe_sra} --stdout "' + f +
                      '" >> {params.sampleID}/{params.sampleID}.fastq.gz')
            else:
                shell('cat "{path_in}/' + f +
                      '" >> {params.sampleID}/{params.sampleID}.fastq.gz')
        shell("touch {output}")
