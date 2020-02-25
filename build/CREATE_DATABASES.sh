!/bin/bash

CORES=8
#BASE=~/WORK/exceRpt/build
BASE=/tmp/exceRpt
PATH_TMP=$BASE/tmp

##
## INSTALL DEPENDENCIES
##
PATH_BIN=$BASE/bin

## run the script to install the dependencies:
./INSTALL_DEPENDENCIES.sh


EXE_EXCERPT_TOOLS=$PATH_BIN/exceRpt/exceRpt_Tools.jar
EXE_THUNDER=$PATH_BIN/Thunder/Thunder.jar
EXE_STAR=$PATH_BIN/STAR/bin/MacOSX_x86_64/STAR
EXE_BOWTIE1=$PATH_BIN/bowtie1/bowtie
EXE_BOWTIE1_BUILD=$PATH_BIN/bowtie1/bowtie-build
EXE_BOWTIE2=$PATH_BIN/bowtie2/bowtie2
EXE_BOWTIE2_BUILD=$PATH_BIN/bowtie2/bowtie2-build
EXE_SAMTOOLS=$PATH_BIN/samtools/samtools
PATH_FASTQC=$PATH_BIN/FastQC
PATH_BBMAP=$PATH_BIN/bbmap
EXE_SRA_TOOLKIT=$PATH_BIN/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump




##
## BUILD DATABASES
##

PATH_DB=$BASE/DATABASE
mkdir -p $PATH_DB
PATH_FA=$BASE/fastaFiles
mkdir -p $PATH_FA



##
## Sync fasta files from S3
##  - these are hard/annoying to find/construct/reformat from online databases
##
aws s3 sync s3://kitchen-mgh-public/exceRpt/DATABASE/fasta_static $PATH_FA






##
## Make pseudo-random data for seeding
##  - this is what was used to create the randomBits.dat seed file
##
#get_seeded_random()
#{
#  seed="$1"
#  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
#    </dev/zero 2>/dev/null
#}
#get_seeded_random 7 > $PATH_DB/randomBits.dat



cd $PATH_FA


##
## Build the STAR index for UniVec and rRNA contaminants
##
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
gzip -c $PATH_FA/UniVec_Core > $PATH_FA/UniVec_Core.contaminants.fa.gz

## hg38
mkdir -p $PATH_DB/hg38/STAR_INDEX_Univec_rRNA
gunzip -c $PATH_FA/hg38_rRNA_45S_5S.fa.gz | sed 's/^>/>rRNA:/' > $PATH_FA/tmp.fa
gunzip -c $PATH_FA/UniVec_Core.contaminants.fa.gz | sed 's/^>/>UniVec:/' >> $PATH_FA/tmp.fa
#cat $PATH_FA/tmp.fa | grep -c "^>" ; cat $PATH_FA/tmp.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg38/STAR_INDEX_Univec_rRNA --genomeFastaFiles $PATH_FA/tmp.fa --genomeSAindexNbases 9 --genomeChrBinNbits 8
rm $PATH_FA/tmp.fa
## hg19
mkdir -p $PATH_DB/hg19/STAR_INDEX_Univec_rRNA
gunzip -c $PATH_FA/hg19_rRNA_45S_5S_extendedregions.fa.gz | sed 's/^>/>rRNA:/' > $PATH_FA/tmp.fa
gunzip -c $PATH_FA/UniVec_Core.contaminants.fa.gz | sed 's/^>/>UniVec:/' >> $PATH_FA/tmp.fa
#cat $PATH_FA/tmp.fa | grep -c "^>" ; cat $PATH_FA/tmp.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg19/STAR_INDEX_Univec_rRNA --genomeFastaFiles $PATH_FA/tmp.fa --genomeSAindexNbases 9 --genomeChrBinNbits 8
rm $PATH_FA/tmp.fa
## mm10
mkdir -p $PATH_DB/mm10/STAR_INDEX_Univec_rRNA
gunzip -c $PATH_FA/mm10_rRNA_all.fa.gz | sed 's/^>/>rRNA:/' > $PATH_FA/tmp.fa
gunzip -c $PATH_FA/UniVec_Core.contaminants.fa.gz | sed 's/^>/>UniVec:/' >> $PATH_FA/tmp.fa
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/mm10/STAR_INDEX_Univec_rRNA --genomeFastaFiles $PATH_FA/tmp.fa --genomeSAindexNbases 9 --genomeChrBinNbits 8
rm $PATH_FA/tmp.fa
rm Log.out



##
## Build the STAR index of the genome
##
## hg38
mkdir -p $PATH_DB/hg38/STAR_INDEX_genome
gunzip -c $PATH_FA/hg38.fa.gz > $PATH_FA/hg38.fa
#$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg38/STAR_INDEX_genome --genomeFastaFiles $PATH_FA/hg38.fa --genomeSAindexNbases 14 --genomeChrBinNbits 18
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg38/STAR_INDEX_genome --genomeFastaFiles $PATH_FA/hg38.fa --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 4
rm $PATH_FA/hg38.fa
## hg19
mkdir -p $PATH_DB/hg19/STAR_INDEX_genome
gunzip -c $PATH_FA/hg19.fa.gz > $PATH_FA/hg19.fa
#$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg19/STAR_INDEX_genome --genomeFastaFiles $PATH_FA/hg19.fa --genomeSAindexNbases 14 --genomeChrBinNbits 18
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg19/STAR_INDEX_genome --genomeFastaFiles $PATH_FA/hg19.fa --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 4
rm $PATH_FA/hg19.fa
## mm10
mkdir -p $PATH_DB/mm10/STAR_INDEX_genome
gunzip -c $PATH_FA/mm10.fa.gz > $PATH_FA/mm10.fa
#$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/mm10/STAR_INDEX_genome --genomeFastaFiles $PATH_FA/mm10.fa --genomeSAindexNbases 14 --genomeChrBinNbits 18
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/mm10/STAR_INDEX_genome --genomeFastaFiles $PATH_FA/mm10.fa --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 4
rm $PATH_FA/mm10.fa



##
## Download and fix miRNA references
##
mkdir -p $PATH_FA/miRNA
cd $PATH_FA/miRNA


VER_MIRBASE=22
wget ftp://mirbase.org/pub/mirbase/22.1/hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/22.1/mature.fa.gz
gunzip -c hairpin.fa.gz > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin.fa
gunzip -c mature.fa.gz > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature.fa

U2T="java -jar $EXE_EXCERPT_TOOLS Fasta_U2T"
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin.fa | sed 's/ /:/g' > $PATH_FA/miRNA/tmp.fa; $U2T $PATH_FA/miRNA/tmp.fa > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_ALL.fa
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature.fa | sed 's/ /:/g' > $PATH_FA/miRNA/tmp.fa; $U2T $PATH_FA/miRNA/tmp.fa > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_ALL.fa
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_ALL.fa | sed 's/U/T/g' | tr '\n' '=' | tr '>' '\n' | grep "^hsa-" | tr '\n' '>' | tr '=' '\n' | sed '1s/^/>/' | grep -v "^>$" > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa.fa
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_ALL.fa | sed 's/U/T/g' | tr '\n' '=' | tr '>' '\n' | grep "^hsa-" | tr '\n' '>' | tr '=' '\n' | sed '1s/^/>/' | grep -v "^>$" > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_hsa.fa
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_ALL.fa | sed 's/U/T/g' | tr '\n' '=' | tr '>' '\n' | grep "^mmu-" | tr '\n' '>' | tr '=' '\n' | sed '1s/^/>/' | grep -v "^>$" > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_mmu.fa
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_ALL.fa | sed 's/U/T/g' | tr '\n' '=' | tr '>' '\n' | grep "^mmu-" | tr '\n' '>' | tr '=' '\n' | sed '1s/^/>/' | grep -v "^>$" > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_mmu.fa

cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_ALL.fa | tr '\n' '=' | tr '>' '\n' | grep -v "^hsa-" | tr '\n' '>' | tr '=' '\n' | sed '$ s/.$//' > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_NONhsa.fa
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_ALL.fa | tr '\n' '=' | tr '>' '\n' | grep -v "^mmu-" | tr '\n' '>' | tr '=' '\n' | sed '$ s/.$//' > $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_NONmmu.fa

rm $PATH_FA/miRNA/tmp.fa


##
## Map hairpin sequences to the respective genomes and then map mature miRNAs to hairpins
##

## Build bowtie2 indices of the genomes for aligning hairpin sequences to
$EXE_BOWTIE2_BUILD -f --threads $CORES $PATH_FA/hg38.fa.gz $PATH_TMP/bowtie2_hg38
$EXE_BOWTIE2_BUILD -f --threads $CORES $PATH_FA/hg19.fa.gz $PATH_TMP/bowtie2_hg19
$EXE_BOWTIE2_BUILD -f --threads $CORES $PATH_FA/mm10.fa.gz $PATH_TMP/bowtie2_mm10

## Map PRECURSOR miRNAs to the genome (bowtie 2 used simply to handle long ‘reads’, turn off indels & splicing)
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa.fa \
	| $EXE_BOWTIE2 -p $CORES -N 0 --rdg 1000,1000 --rfg 1000,1000 --all --reorder --no-head --end-to-end -f \
				  	-x $PATH_TMP/bowtie2_hg38 -U - > $PATH_DB/hg38/miRNA_precursor2genome.sam
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa.fa \
        | $EXE_BOWTIE2 -p $CORES -N 0 --rdg 1000,1000 --rfg 1000,1000 --all --reorder --no-head --end-to-end -f \
                                        -x $PATH_TMP/bowtie2_hg19 -U - > $PATH_DB/hg19/miRNA_precursor2genome.sam
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_mmu.fa \
        | $EXE_BOWTIE2 -p $CORES -N 0 --rdg 1000,1000 --rfg 1000,1000 --all --reorder --no-head --end-to-end -f \
                                        -x $PATH_TMP/bowtie2_mm10 -U - > $PATH_DB/mm10/miRNA_precursor2genome.sam

## Create bowtie1 references for the hairpin sequences
$EXE_BOWTIE1_BUILD $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa.fa $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa
$EXE_BOWTIE1_BUILD $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_mmu.fa $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_mmu

## Map MATURE miRNAs to the HAIRPINs (bowtie 1 ungapped)
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_hsa.fa \
	| $EXE_BOWTIE1 -p $CORES -v 1 --all --sam --sam-nohead --best --strata \
	-f $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa - > $PATH_DB/hg38/miRNA_mature2precursor.sam
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_hsa.fa \
        | $EXE_BOWTIE1 -p $CORES -v 1 --all --sam --sam-nohead --best --strata \
        -f $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa - > $PATH_DB/hg19/miRNA_mature2precursor.sam
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_mature_mmu.fa \
        | $EXE_BOWTIE1 -p $CORES -v 1 --all --sam --sam-nohead --best --strata \
        -f $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_mmu - > $PATH_DB/mm10/miRNA_mature2precursor.sam







##
## Download and fix tRNA references
##
cd $PATH_FA

## hg38
wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa
cat $PATH_FA/hg38-tRNAs.fa | sed 's/Homo_sapiens_//g' | awk -F " " '{print $1}' > $PATH_FA/hg38-tRNAs_modifiedHeaders.fa
## hg19
wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.fa
cat $PATH_FA/hg19-tRNAs.fa | sed 's/Homo_sapiens_//g' | awk -F " " '{print $1}' > $PATH_FA/hg19-tRNAs_modifiedHeaders.fa
## mm10
wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.fa
cat $PATH_FA/mm10-tRNAs.fa | sed 's/Mus_musculus_//g' | awk -F " " '{print $1}' > $PATH_FA/mm10-tRNAs_modifiedHeaders.fa

rm $PATH_FA/hg38-tRNAs.fa
rm $PATH_FA/hg19-tRNAs.fa
rm $PATH_FA/mm10-tRNAs.fa

gzip $PATH_FA/hg38-tRNAs_modifiedHeaders.fa
gzip $PATH_FA/hg19-tRNAs_modifiedHeaders.fa
gzip $PATH_FA/mm10-tRNAs_modifiedHeaders.fa 




##
## Download and fix gencode references
##
cd $PATH_FA

VER_GENCODE_HSA=33
VER_GENCODE_MMU=M24

GTF2FA="java -Xmx20G -jar $EXE_THUNDER GTF2Fasta"

## hg38
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$VER_GENCODE_HSA/gencode.v$VER_GENCODE_HSA.annotation.gtf.gz
gunzip -c $PATH_FA/gencode.v$VER_GENCODE_HSA.annotation.gtf.gz > $PATH_DB/hg38/gencodeAnnotation.gtf
cp $PATH_FA/gencode.v$VER_GENCODE_HSA.annotation.gtf.gz $PATH_DB/hg38/gencodeAnnotation.gtf.gz
gunzip -c $PATH_FA/hg38.fa.gz > $PATH_FA/hg38.fa
$GTF2FA -a $PATH_DB/hg38/gencodeAnnotation.gtf -i $PATH_FA/hg38.fa -N false -o $PATH_FA/hg38_gencode.v$VER_GENCODE_HSA.fa -t transcript_type,transcript_name
gzip $PATH_FA/hg38_gencode.v$VER_GENCODE_HSA.fa
rm $PATH_DB/hg38/gencodeAnnotation.gtf

## hg19
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip -c $PATH_FA/gencode.v19.annotation.gtf.gz > $PATH_DB/hg19/gencodeAnnotation.gtf
cp $PATH_FA/gencode.v19.annotation.gtf.gz $PATH_DB/hg19/gencodeAnnotation.gtf.gz
gunzip -c $PATH_FA/hg19.fa.gz > $PATH_FA/hg19.fa
$GTF2FA -a $PATH_DB/hg38/gencodeAnnotation.gtf -i $PATH_FA/hg19.fa -N false -o $PATH_FA/hg19_gencode.v19.fa -t transcript_type,transcript_name
gzip $PATH_FA/hg19_gencode.v19.fa
rm $PATH_DB/hg19/gencodeAnnotation.gtf

## mm10
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_$VER_GENCODE_MMU/gencode.v$VER_GENCODE_MMU.annotation.gtf.gz
gunzip -c gencode.v$VER_GENCODE_MMU.annotation.gtf.gz > $PATH_DB/mm10/gencodeAnnotation.gtf
cp $PATH_FA/gencode.v$VER_GENCODE_MMU.annotation.gtf.gz $PATH_DB/mm10/gencodeAnnotation.gtf.gz
gunzip -c $PATH_FA/mm10.fa.gz > $PATH_FA/mm10.fa
$GTF2FA -a $PATH_DB/mm10/gencodeAnnotation.gtf -i $PATH_FA/mm10.fa -N false -o $PATH_FA/mm10_gencode.v$VER_GENCODE_MMU.fa -t transcript_type,transcript_name
gzip $PATH_FA/mm10_gencode.v$VER_GENCODE_MMU.fa
rm $PATH_DB/mm10/gencodeAnnotation.gtf




##
## Download and fix piRNA references
##

##
## piRNAbank - THIS IS NOW DEFUNCT!
##
## hg38
#wget http://pirnabank.ibab.ac.in/Human.tar.gz
#tar -xf Human.tar.gz
#cat Human/hsa_piR_*.txt > piRNA_piRNABank_human.fa
#java -Xmx10G -jar ~/bin/smallRNAPipeline/exceRpt_Tools.jar RemoveFastaDuplicates -o piRNA_piRNABank_human_noDups.fa -s piRNA_piRNABank_human.fa
#rm -r Human

## mm10
#wget http://pirnabank.ibab.ac.in/Mouse.tar.gz
#tar -xf Mouse.tar.gz
#cat Mouse/mmu_piR_*.txt > piRNA_piRNABank_mouse.fa
#java -Xmx10G -jar ~/bin/smallRNAPipeline/exceRpt_Tools.jar RemoveFastaDuplicates -o piRNA_piRNABank_mouse_noDups.fa -s piRNA_piRNABank_mouse.fa
#rm -r Mouse

##
## piRBase
##
## --> NOT going to do anything further with this as there are WAY too many sequences here - can't possibly be real
##
#mkdir -p $PATH_FA/piRNA/piRBase
#cd $PATH_FA/piRNA/piRBase
#wget http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/fasta/hsa.fa.gz
#wget http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/fasta/mmu.fa.gz
#java -Xmx10G -jar ~/bin/smallRNAPipeline/exceRpt_Tools.jar RemoveFastaDuplicates -o piRNA_piRNABank_mouse_noDups.fa -s piRNA_piRNABank_mouse.fa


##
## NEW piRNA reference
## 
## HUMAN
## - https://rnacentral.org/search?q=piRNA*%20AND%20TAXONOMY:%229606%22%20AND%20expert_db:%22ENA%22%20AND%20rna_type:%22piRNA%22%20AND%20has_genomic_coordinates:%22True%22
##
## MOUSE
## - https://rnacentral.org/search?q=piRNA*%20AND%20rna_type:%22piRNA%22%20AND%20has_genomic_coordinates:%22True%22%20AND%20TAXONOMY:%2210090%22%20AND%20expert_db:%22ENA%22%20AND%20length:%5B17%20TO%2050%5D
##
mkdir -p $PATH_FA/piRNA/RNAcentral
cd $PATH_FA/piRNA/RNAcentral

RMDUP="java -Xmx10G -jar $EXE_EXCERPT_TOOLS RemoveFastaDuplicates"

gunzip -c $PATH_FA/piRNA/RNAcentral/piRNA_AND_TAXONOMY9606_AND_expert_dbENA_AND_rna_typepiRNA_AND_has_genomic_coordinatesTrue.fasta.gz > $PATH_FA/piRNA/RNAcentral/tmp.fa
$RMDUP -o $PATH_FA/piRNA/RNAcentral/piRNA_RNAcentral_human.fa -s $PATH_FA/piRNA/RNAcentral/tmp.fa
gzip $PATH_FA/piRNA/RNAcentral/piRNA_RNAcentral_human.fa
rm $PATH_FA/piRNA/RNAcentral/tmp.fa

gunzip -c $PATH_FA/piRNA/RNAcentral/piRNA_AND_rna_typepiRNA_AND_has_genomic_coordinatesTrue_AND_TAXONOMY10090_AND_expert_dbENA_AND_length17_TO_50.fasta.gz > $PATH_FA/piRNA/RNAcentral/tmp.fa
$RMDUP -o $PATH_FA/piRNA/RNAcentral/piRNA_RNAcentral_mouse.fa -s $PATH_FA/piRNA/RNAcentral/tmp.fa
gzip $PATH_FA/piRNA/RNAcentral/piRNA_RNAcentral_mouse.fa
rm $PATH_FA/piRNA/RNAcentral/tmp.fa



##
## Circular RNAs from CircBase
##
mkdir -p $PATH_FA/circRNA/circBase
cd $PATH_FA/circRNA/circBase

wget http://www.circbase.org/download/human_hg19_circRNAs_putative_spliced_sequence.fa.gz
wget http://www.circbase.org/download/mouse_mm9_circRNAs_putative_spliced_sequence.fa.gz
#wget http://www.circbase.org/download/hsa_hg19_circRNA.txt
#wget http://www.circbase.org/download/mmu_mm9_circRNA.txt



##
## Build the STAR index of all libraries
##
## hg38
mkdir -p $PATH_DB/hg38/STAR_INDEX_transcriptome
## concatenate the individual libraries into a single fasta
cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_hsa.fa | sed 's/^>/>miRNA:/' > $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/hg38-tRNAs_modifiedHeaders.fa.gz | sed 's/^>/>tRNA:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/piRNA/RNAcentral/piRNA_RNAcentral_human.fa.gz | sed 's/^>/>piRNA:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/hg38_gencode.v$VER_GENCODE_HSA.fa.gz | sed 's/^>/>gencode:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/circRNA/circBase/human_hg19_circRNAs_putative_spliced_sequence.fa.gz | sed 's/^>/>circRNA:/' >> $PATH_FA/tmp.allLibs.fa
## Count the number of references and the number of nucleotides
#cat $PATH_FA/tmp.allLibs.fa | grep -c "^>" ; cat $PATH_FA/tmp.allLibs.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg38/STAR_INDEX_transcriptome --genomeFastaFiles $PATH_FA/tmp.allLibs.fa --genomeSAindexNbases 14 --genomeChrBinNbits 11

## hg19
mkdir -p $PATH_DB/hg19/STAR_INDEX_transcriptome
## concatenate the individual libraries into a single fasta
cat $PATH_FA/miRBase_v$VER_MIRBASE\_hairpin_hsa.fa | sed 's/^>/>miRNA:/' > $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/hg19-tRNAs_modifiedHeaders.fa.gz | sed 's/^>/>tRNA:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/piRNA/RNAcentral/piRNA_RNAcentral_human.fa.gz | sed 's/^>/>piRNA:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/hg19_gencode.v19.fa.gz | sed 's/^>/>gencode:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/circRNA/circBase/human_hg19_circRNAs_putative_spliced_sequence.fa.gz | sed 's/^>/>circRNA:/' >> $PATH_FA/tmp.allLibs.fa
## Count the number of references and the number of nucleotides
#cat $PATH_FA/tmp.allLibs.fa | grep -c "^>" ; cat $PATH_FA/tmp.allLibs.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg19/STAR_INDEX_transcriptome --genomeFastaFiles $PATH_FA/tmp.allLibs.fa --genomeSAindexNbases 14 --genomeChrBinNbits 11

## mm10
mkdir -p $PATH_DB/mm10/STAR_INDEX_transcriptome
## concatenate the individual libraries into a single fasta
cat $PATH_FA/miRBase_v$VER_MIRBASE\_hairpin_mmu.fa | sed 's/^>/>miRNA:/' > $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/mm10-tRNAs_modifiedHeaders.fa.gz | sed 's/^>/>tRNA:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/piRNA/RNAcentral/piRNA_RNAcentral_mouse.fa.gz | sed 's/^>/>piRNA:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/mm10_gencode.v$VER_GENCODE_MMU.fa.gz | sed 's/^>/>gencode:/' >> $PATH_FA/tmp.allLibs.fa
gunzip -c $PATH_FA/circRNA/circBase/mouse_mm9_circRNAs_putative_spliced_sequence.fa.gz | sed 's/^>/>circRNA:/' >> $PATH_FA/tmp.allLibs.fa
## Count the number of references and the number of nucleotides
#cat $PATH_FA/tmp.allLibs.fa | grep -c "^>" ; cat $PATH_FA/tmp.allLibs.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/mm10/STAR_INDEX_transcriptome --genomeFastaFiles $PATH_FA/tmp.allLibs.fa --genomeSAindexNbases 14 --genomeChrBinNbits 11

rm $PATH_FA/tmp.allLibs.fa



##
## Build the STAR index for the repetitive elements
##
## hg38
mkdir -p $PATH_DB/hg38/STAR_INDEX_repetitiveElements
gunzip -c $PATH_FA/hg38.RE.fixedHeaders.fa.gz > $PATH_FA/hg38.RE.fixedHeaders.fa
#cat $FASTA/hg38.RE.fixedHeaders.fa | grep -c "^>" ; cat $FASTA/hg38.RE.fixedHeaders.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg38/STAR_INDEX_repetitiveElements --genomeFastaFiles $PATH_FA/hg38.RE.fixedHeaders.fa --genomeSAindexNbases 14 --genomeChrBinNbits 8 --genomeSAsparseD 4
rm $PATH_FA/hg38.RE.fixedHeaders.fa
## hg19
mkdir -p $PATH_DB/hg19/STAR_INDEX_repetitiveElements
gunzip -c $PATH_FA/hg19.RE.fixedHeaders.fa.gz > $PATH_FA/hg19.RE.fixedHeaders.fa
#cat $FASTA/hg19.RE.fixedHeaders.fa | grep -c "^>" ; cat $FASTA/hg19.RE.fixedHeaders.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/hg19/STAR_INDEX_repetitiveElements --genomeFastaFiles $PATH_FA/hg19.RE.fixedHeaders.fa --genomeSAindexNbases 14 --genomeChrBinNbits 8 --genomeSAsparseD 4
rm $PATH_FA/hg19.RE.fixedHeaders.fa
## mm10
mkdir -p $PATH_DB/mm10/STAR_INDEX_repetitiveElements
gunzip -c $PATH_FA/mm10.RE.fixedHeaders.fa.gz > $PATH_FA/mm10.RE.fixedHeaders.fa
#cat $FASTA/mm10.RE.fixedHeaders.fa | grep -c "^>" ; cat $FASTA/mm10.RE.fixedHeaders.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/mm10/STAR_INDEX_repetitiveElements --genomeFastaFiles $PATH_FA/mm10.RE.fixedHeaders.fa --genomeSAindexNbases 14 --genomeChrBinNbits 8 --genomeSAsparseD 4
rm $PATH_FA/mm10.RE.fixedHeaders.fa





##
## Build the STAR index of all miRBase miRNAs
##
mkdir -p $PATH_DB/miRBase/STAR_INDEX_miRBaseAll
#cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_ALL.fa | grep -c "^>"; cat $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_ALL.fa | grep -v "^>" | wc -c
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/miRBase/STAR_INDEX_miRBaseAll --genomeFastaFiles $PATH_FA/miRNA/miRBase_v$VER_MIRBASE\_hairpin_ALL.fa --genomeSAindexNbases 10 --genomeChrBinNbits 7



##
## Build the STAR index of the ribosomal RNAs
##
VERSION_RDP = 11

mkdir -p $PATH_FA/ribosomeDB
cd $PATH_FA/ribosomeDB
wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Archaea_unaligned.fa.gz
wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Fungi_unaligned.fa.gz
gunzip -c $PATH_FA/ribosomeDB/current_Archaea_unaligned.fa.gz | tr '[:blank:]' '_' > $PATH_FA/ribosomeDB/current_ribosomeDB_unaligned.fa
gunzip -c $PATH_FA/ribosomeDB/current_Bacteria_unaligned.fa.gz | tr '[:blank:]' '_' >> $PATH_FA/ribosomeDB/current_ribosomeDB_unaligned.fa
gunzip -c $PATH_FA/ribosomeDB/current_Fungi_unaligned.fa.gz | tr '[:blank:]' '_' >> $PATH_FA/ribosomeDB/current_ribosomeDB_unaligned.fa
$EXE_STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $PATH_DB/ribosomeDatabase/exogenous_rRNAs --genomeFastaFiles $PATH_FA/ribosomeDB/current_ribosomeDB_unaligned.fa --genomeSAindexNbases 10 --genomeChrBinNbits 7
rm $PATH_FA/ribosomeDB/current_ribosomeDB_unaligned.fa





##
## Compress the DB
##

echo -e "$PATH_DB/STAR_Parameters_Endogenous_smallRNA.in"  > $PATH_TMP/filesToCompress_generic.txt
echo -e "$PATH_DB/adapters.fa" 				  >> $PATH_TMP/filesToCompress_generic.txt
echo -e "$PATH_DB/phiX.fa" 				  >> $PATH_TMP/filesToCompress_generic.txt
echo -e "$PATH_DB/randomBits.dat" 			  >> $PATH_TMP/filesToCompress_generic.txt

cp $PATH_TMP/filesToCompress_generic.txt $PATH_TMP/filesToCompress_hg19.txt
echo -e "$PATH_DB/hg19" >> $PATH_TMP/filesToCompress_hg19.txt

cp $PATH_TMP/filesToCompress_generic.txt $PATH_TMP/filesToCompress_hg38.txt
echo -e "$PATH_DB/hg38" >> $PATH_TMP/filesToCompress_hg38.txt

cp $PATH_TMP/filesToCompress_generic.txt $PATH_TMP/filesToCompress_mm10.txt
echo -e "$PATH_DB/mm10" >> $PATH_TMP/filesToCompress_mm10.txt


tar -cv -T $PATH_TMP/filesToCompress_hg19.txt -f $PATH_DB/exceRptDB_v5_hg19.tar
tar -cvf exceRptDB_v5_hg38.tgz hg38
tar -cvf exceRptDB_v5_mm10.tgz mm10
#tar -cvz -T filesToCompress_EXO_miRNArRNA.txt -f exceRptDB_v4_EXOmiRNArRNA.tgz



##
## Sync with S3
##
## Sync fasta files from S3
aws s3 sync $PATH_DB s3://kitchen-mgh-public/exceRpt/DATABASE/v5.0




##
## TEST
##
mkdir -p $BASE/testOutput
make -nf $PATH_BIN/exceRpt/exceRpt_smallRNA \
    INPUT_FILE_PATH=$PATH_BIN/exceRpt/ExampleData/testData_human.fastq.gz \
    OUTPUT_DIR=$BASE/testOutput \
    EXE_DIR=$PATH_BIN \
    DATABASE_PATH=$PATH_DB \
    N_THREADS=16


EXE_DIR=$PATH_BIN




