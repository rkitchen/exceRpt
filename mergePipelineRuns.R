#######################################################################################
##                                                                                   ##
## Script to combine pipeline runs for individual samples into something more useful ##
##                                                                                   ##
## Author: Rob Kitchen (rob.kitchen@yale.edu)                                        ##
##                                                                                   ##
## Version 2.0.1 (2015-01-18)                                                        ##
##                                                                                   ##
#######################################################################################


##
## Define inputs:
##
args<-commandArgs(TRUE)
data.dir = args[1]
output.dir = args[2]


##
## Check inputs
##
## if no data directory is specified, throw an error message and try a hard coded path (for testing)
if(is.na(data.dir)){
  cat("ERROR: no input data directory specified!\n\n")
  cat("Usage: R CMD BATCH mergePipelineRuns_PRODUCTION.R <data path> [output path]\n\n")
  #data.dir = "~/Box Sync/Work for other people/FraminghamSmallRNA/Data"
  #output.dir = "~/Box Sync/Work for other people/FraminghamSmallRNA"
  #data.dir = "/Users/robk/WORK/YALE_offline/miRNA/DataSets"
  #output.dir = "/Users/robk/Box Sync/Work/miRNA/DataSets"
  #data.dir = "/Users/robk/WORK/YALE_offline/exRNA/Galas_testData"
  #data.dir = "/Users/robk/WORK/YALE_offline/exRNA/newPipelineTest"	
  #data.dir = "/Users/robk/WORK/YALE_offline/exRNA/AlCharest/results"	
}

## if not output.dir, default back to data.dir
if(is.na(output.dir)){
  output.dir = data.dir
}



##
## check dependencies
##
if(!"plyr" %in% rownames(installed.packages())) { install.packages("plyr") }
if(!"gplots" %in% rownames(installed.packages())) { install.packages("gplots") }
if(!"marray" %in% rownames(installed.packages())) { source("http://bioconductor.org/biocLite.R"); biocLite("marray") }
if(!"reshape2" %in% rownames(installed.packages())) { install.packages("reshape2") }
if(!"ggplot2" %in% rownames(installed.packages())) { install.packages("ggplot2") }
require(plyr)
require(gplots)
require(marray)
require(reshape2)
require(ggplot2)


######################################################################################
######################################################################################
######################################################################################
##             There should be no need to change the code below...                  ##
######################################################################################
######################################################################################
######################################################################################


##
## Function to recursively search a given directory for pipeline output
##
SearchForSampleData = function(base.dir, directory=""){
  to.return = NULL
  dir.use = paste(base.dir, directory, sep = "/")
  subdirs = dir(dir.use)
  if(length(subdirs) > 0){
    i.stats = grep(".stats$", subdirs, perl=T)
    i.zip = grep(".zip$", subdirs, perl=T)
    
    ## handle decompressed pipeline output
    if(length(i.stats) > 0){
      tmp = gsub(".stats","",subdirs[i.stats])
      to.return = c(to.return, paste(dir.use,tmp,sep="/"))
    }
    
    ## handle zipped pipeline output
    if(length(i.zip) > 0){
      for(x in i.zip){
        tmp.contents = unzip(paste(dir.use,subdirs[x],sep="/"), list=T)[,1]
        if(length(grep(".stats$", tmp.contents, perl=T)) > 0){
          try(unzip(paste(dir.use,subdirs[x],sep="/"), exdir=gsub(".zip","",paste(dir.use,subdirs[x],sep="/")), overwrite=FALSE), silent=T)
          to.return = c(to.return, paste(dir.use,gsub(".zip$","",subdirs[x]),sep="/",gsub(".stats$","",tmp.contents[grep(".stats$", tmp.contents)])))
        }
      }
    }
    
    ## handle unknown directories
    i.known = c(i.stats,i.zip)
    if(length(i.known) == 0){		# there are no .stats or .zip files in the directory!
      i.unknown = 1:length(subdirs)
    }else{
      i.unknown = (1:length(subdirs))[-i.known]
    }
    
    if(length(i.unknown) > 0){
      for(x in subdirs[i.unknown]){
        to.return = c(to.return, SearchForSampleData(dir.use, x))
      }
    }
    
    return(unique(to.return))
  }
}


##
## Get list of paths containing pipeline output
##  -- Kill script if we do not have any samples to process
##
samplePathList = SearchForSampleData(data.dir,"")
NumberOfCompatibleSamples = length(samplePathList)
stopifnot(NumberOfCompatibleSamples > 0)


##
## Function to read the sample data
##
readData = function(path,file,countColumnIndex=4){
  if(file %in% dir(path)){
    tmp = read.table(paste(path,file,sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F, fill=T, colClasses=c("character","numeric","numeric","numeric","numeric","numeric"))
    return(ddply(tmp,"name",function(block){ coords=""; if(ncol(block)==7){ coords=paste(block[,7],collapse="_") }; return(c(sum(block[,countColumnIndex]), coords)) }))
  }else{
    return(NULL)
  }
}


##
## Create objects to contain the data
##
sample.data = vector(mode="list",length=length(samplePathList))
allIDs.miRNA = NULL
allIDs.tRNA = NULL
allIDs.piRNA = NULL
allIDs.gencode = NULL
allIDs.exogenous_miRNA = NULL
mapping.stats = matrix(0,nrow=length(samplePathList),ncol=16, dimnames=list(1:length(samplePathList), c("input","clipped","calibrator","UniVec_contaminants","rRNA","not_rRNA","genome","miRNA_sense","miRNA_antisense","tRNA_sense","tRNA_antisense","piRNA_sense","piRNA_antisense","Gencode_sense","Gencode_antisense","miRNA_exogenous_sense")))
read.lengths = matrix(0,nrow=length(samplePathList),ncol=76,dimnames=list(1:length(samplePathList), 0:75))


##
## Loop through all samples and read the pipeline output
##
for(i in 1:length(samplePathList)){
  ## Parse the sampleID from the path:
  tmp = unlist(strsplit(samplePathList[i], "/"))
  thisSampleID = tmp[length(tmp)]
  
  ## Read sample mapping stats
  tmp.stats = read.table(paste(samplePathList[i],".stats",sep=""), stringsAsFactors=F, fill=T, header=T, sep="\t",skip=0)
  mapping.stats[i, match(tmp.stats[,1], colnames(mapping.stats))] = as.numeric(tmp.stats[,2])
  rownames(mapping.stats)[i] = thisSampleID
  
  ##
  ## Read the clipped read lengths  
  ##
  if(length(grep("clipped.readLengths.txt$", dir(samplePathList[i]))) == 1){
    tmp = read.table(paste(samplePathList[i], dir(samplePathList[i])[grep("clipped.readLengths.txt$", dir(samplePathList[i]))], sep="/"))
    read.lengths[i, 1:ncol(tmp)] = as.numeric(tmp[2,])
    rownames(read.lengths)[i] = thisSampleID
  }
  ## TODO:  Read the read-lengths input to sRNAbench in /stat?
  if("readLengthAnalysis.txt" %in% dir(paste(samplePathList[i],"/stat",sep=""))){
    #read.table(paste(samplePathList[i],"/stat/readLengthAnalysis.txt",sep=""), sep="\t", header=T)
  }
  
  
  ##
  ## Read the adapter sequence
  ##
  if(paste(thisSampleID,".adapterSeq",sep="") %in% dir(samplePathList[i])){
    adapterSeq = as.character(read.table(paste(samplePathList[i],"/",thisSampleID,".adapterSeq",sep=""))[1,1])
  }
  
  
  ##
  ## Read sample data
  ##
  if("mature_sense_singleA.grouped" %in% dir(samplePathList[i])){
    ## For the NEW version of sRNABench:
    countColumnIndex=4
    miRNA = readData(samplePathList[i],"mature_sense.grouped", countColumnIndex)
    tRNA = readData(samplePathList[i],"mm10_gencodeM4_tRNAs_sense.grouped", countColumnIndex)
    piRNA = readData(samplePathList[i],"mm10_piRNAs_sense.grouped", countColumnIndex)
    gencode = readData(samplePathList[i],"mm10_gencodeM4_sense.grouped", countColumnIndex)
  }else{
    ## Not valid pipeline output!
  }
  
  # Update list of detected smallRNA IDs
  allIDs.miRNA = unique(c(allIDs.miRNA, as.character(miRNA[,1])))
  allIDs.tRNA = unique(c(allIDs.tRNA, as.character(tRNA[,1])))
  allIDs.piRNA = unique(c(allIDs.piRNA, as.character(piRNA[,1])))
  allIDs.gencode = unique(c(allIDs.gencode, as.character(gencode[,1])))
  
  sample.data[[i]] = list("miRNA"=miRNA, "tRNA"=tRNA, "piRNA"=piRNA, "gencode"=gencode, adapterSeq=adapterSeq)
  names(sample.data)[i] = thisSampleID
  
  
  ##
  ## Read plant/virus alignments (if applicable)
  ##
  if("EXOGENOUS_miRNA" %in% dir(samplePathList[i])){
    tmp.dir = paste(samplePathList[i],"EXOGENOUS_miRNA",sep="/")
	if("mature_sense_singleA.grouped" %in% dir(tmp.dir)){
      ## For the NEW version of sRNABench:
      sample.data[[i]]$exogenous_miRNA = readData(tmp.dir,"mature_sense.grouped", 4)
	  allIDs.exogenous_miRNA = unique(c(allIDs.exogenous_miRNA, as.character(sample.data[[i]]$exogenous_miRNA[,1])))
    }
  }
}


##
## Save result
##
allIDs = list("miRNA"=allIDs.miRNA, "tRNA"=allIDs.tRNA, "piRNA"=allIDs.piRNA, "gencode"=allIDs.gencode, "exogenous_miRNA"=allIDs.exogenous_miRNA)
#save(sample.data, mapping.stats, allIDs, file=paste(output.dir,"RawFileData.RData",sep="/"))


##
## Convert sample data to large per-smallRNA expression matrices
##
exprs.miRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$miRNA), dimnames=list(allIDs$miRNA, names(sample.data)))
exprs.tRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$tRNA), dimnames=list(allIDs$tRNA, names(sample.data)))
exprs.piRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$piRNA), dimnames=list(allIDs$piRNA, names(sample.data)))
exprs.gencode = matrix(0,ncol=length(sample.data),nrow=length(allIDs$gencode), dimnames=list(allIDs$gencode, names(sample.data)))
exprs.exogenous_miRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$exogenous_miRNA), dimnames=list(allIDs$exogenous_miRNA, names(sample.data)))
for(i in 1:length(sample.data)){
  exprs.miRNA[match(sample.data[[i]]$miRNA[,1], rownames(exprs.miRNA)),i] = as.numeric(sample.data[[i]]$miRNA[,2])
  exprs.tRNA[match(sample.data[[i]]$tRNA[,1], rownames(exprs.tRNA)),i] = as.numeric(sample.data[[i]]$tRNA[,2])
  exprs.piRNA[match(sample.data[[i]]$piRNA[,1], rownames(exprs.piRNA)),i] = as.numeric(sample.data[[i]]$piRNA[,2])
  exprs.gencode[match(sample.data[[i]]$gencode[,1], rownames(exprs.gencode)),i] = as.numeric(sample.data[[i]]$gencode[,2])
  exprs.exogenous_miRNA[match(sample.data[[i]]$exogenous_miRNA[,1], rownames(exprs.exogenous_miRNA)),i] = as.numeric(sample.data[[i]]$exogenous_miRNA[,2])
}


##
## Remove samples with no miRNA counts
##
#failedSamples = which(colSums(exprs.miRNA) == 0)
#if(length(failedSamples) > 0){
#  #sample.decode[failedSamples, ]
#  #sample.decode = sample.decode[-failedSamples, ]
#  exprs.miRNA = exprs.miRNA[, -failedSamples]
#  exprs.tRNA = exprs.tRNA[, -failedSamples]
#  exprs.piRNA = exprs.piRNA[, -failedSamples]
#  exprs.snoRNA = exprs.snoRNA[, -failedSamples]
#  exprs.Rfam = exprs.Rfam[, -failedSamples]
#  exprs.plantVirus = exprs.plantVirus[, -failedSamples]
#}
dim(exprs.miRNA)
dim(exprs.tRNA)
dim(exprs.piRNA)
dim(exprs.gencode)
dim(exprs.exogenous_miRNA)


##
## Calculate the total number of mapped reads to the rRNA, genome, and exogenous sequences
##
mapping.stats[is.na(mapping.stats)] = 0
mapping.stats = as.data.frame(mapping.stats)
libSizes = list()
libSizes$clipped = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("clipped")])
libSizes$all = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("rRNA","genome","miRNA_exogenous_sense")])
libSizes$endogenous = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("rRNA","genome")])
libSizes$genome = mapping.stats[,colnames(mapping.stats) %in% "genome"]
libSizes$smRNA = mapping.stats[,grep("sense",colnames(mapping.stats))]
libSizes$miRNA = colSums(exprs.miRNA)


##
## Save the raw count data
##
save(exprs.miRNA, exprs.tRNA, exprs.piRNA, exprs.gencode, exprs.exogenous_miRNA, mapping.stats, libSizes, read.lengths, file=paste(output.dir, "smallRNAQuants_ReadCounts.RData", sep="/"))
write.table(exprs.miRNA, file=paste(output.dir, "miRNA_Quantifications_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
write.table(read.lengths, file=paste(output.dir, "ReadLengths.txt", sep="/"), sep="\t", col.names=NA, quote=F)


##
## Calculate the fractions of reads mapping at each stage
##
write.table(mapping.stats, file=paste(output.dir,"readMappingSummary.txt",sep="/"), sep="\t", col.names=NA, quote=F)

fractions = NULL
quals = t(mapping.stats)

if("input" %in% colnames(mapping.stats)){
  inputReads = mapping.stats[,colnames(mapping.stats) %in% c("input")]
}else{
  inputReads = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("rRNA","genome")])
}
if("clipped" %in% colnames(mapping.stats)){
  fractions$clipped = 100*as.numeric(mapping.stats[,colnames(mapping.stats) %in% "clipped"]) / inputReads
}
if("calibrator" %in% colnames(mapping.stats)){
  if( ! is.na(as.numeric(mapping.stats[,colnames(mapping.stats) %in% "calibrator"])[1])){
    fractions$calibrator = 100*as.numeric(mapping.stats[,grep("calibrator",colnames(mapping.stats))]) / inputReads
  }
}
if("UniVec_contaminants" %in% colnames(mapping.stats)){
  fractions$UniVec_contaminants = 100*as.numeric(mapping.stats[,colnames(mapping.stats) %in% "UniVec_contaminants"]) / inputReads
}
if("rRNA" %in% colnames(mapping.stats)){
  fractions$rRNA = 100*as.numeric(mapping.stats[,colnames(mapping.stats) %in% "rRNA"]) / inputReads
}
if("genome" %in% colnames(mapping.stats)){
  fractions$genome = 100*as.numeric(mapping.stats[,colnames(mapping.stats) %in% "genome"]) / inputReads
}
if("miRNA_sense" %in% colnames(mapping.stats)){
  miRNA.counts = mapping.stats[,colnames(mapping.stats) %in% "miRNA_sense"]
  if("miRNA_antisense" %in% colnames(mapping.stats)){
    miRNA.counts = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("miRNA_sense","miRNA_antisense")])
  }
  fractions$miRNA = 100 * miRNA.counts / mapping.stats[,colnames(mapping.stats) %in% "genome"]
  
  mapping.stats = cbind(mapping.stats, miRNA=miRNA.counts)
}

## Finally, convert the computed alignment fractions to a data.frame
fracs = t(as.data.frame(fractions))
colnames(fracs) = rownames(mapping.stats)


##
## Determine heatmap quality colours (green, orange, red) based on read mapping stats
##
quals = fracs

## For the fraction of successfully clipped reads
if("clipped" %in% names(fractions)){
  i.good = which(fractions$clipped > 90)
  i.warn = which(fractions$clipped <= 90  &  fractions$clipped > 50)
  i.bad = which(fractions$clipped <= 50)
  tmp.row = rownames(quals) %in% "clipped"
  quals[tmp.row,i.good] = 100; quals[tmp.row,i.warn] = 50; quals[tmp.row,i.bad] = 0; 
}

## For the calibrator oligo signal
if("calibrator" %in% names(fractions)){
  i.good = which(fractions$calibrator < 10)
  i.warn = which(fractions$calibrator >= 10  &  fractions$calibrator < 50)
  i.bad = which(fractions$calibrator >= 50)
  tmp.row=rownames(quals) %in% "calibrator"
  quals[tmp.row,i.good] = 100; quals[tmp.row,i.warn] = 50; quals[tmp.row,i.bad] = 0; 
}

## For the UniVec contamination
if("UniVec_contaminants" %in% names(fractions)){
  #pass=10; fail=50;
  pass=25; fail=75;
  i.good = which(fractions$UniVec_contaminants < pass)
  i.warn = which(fractions$UniVec_contaminants >= pass  &  fractions$UniVec_contaminants < fail)
  i.bad = which(fractions$UniVec_contaminants >= fail)
  tmp.row=rownames(quals) %in% "UniVec_contaminants"
  quals[tmp.row,i.good] = 100; quals[tmp.row,i.warn] = 50; quals[tmp.row,i.bad] = 0; 
}

## For the rRNA contamination
if("rRNA" %in% names(fractions)){
  #pass=10; fail=50;
  pass=25; fail=75;
  i.good = which(fractions$rRNA < pass)
  i.warn = which(fractions$rRNA >= pass  &  fractions$rRNA < fail)
  i.bad = which(fractions$rRNA >= fail)
  tmp.row=rownames(quals) %in% "rRNA"
  quals[tmp.row,i.good] = 100; quals[tmp.row,i.warn] = 50; quals[tmp.row,i.bad] = 0; 
}

## For genome mapping rate
if("genome" %in% names(fractions)){
  pass=50; fail=25;
  i.good = which(fractions$genome > pass)
  i.warn = which(fractions$genome <= pass  &  fractions$genome > fail)
  i.bad = which(fractions$genome <= fail)
  tmp.row=rownames(quals) %in% "genome"
  quals[tmp.row,i.good] = 100; quals[tmp.row,i.warn] = 50; quals[tmp.row,i.bad] = 0; 
}

## For miRNA mapping rate
if("miRNA" %in% names(fractions)){
  pass=50; fail=25;
  i.good = which(fractions$miRNA > pass)
  i.warn = which(fractions$miRNA <= pass  &  fractions$miRNA > fail)
  i.bad = which(fractions$miRNA <= fail)
  tmp.row=rownames(quals) %in% "miRNA"
  quals[tmp.row,i.good] = 100; quals[tmp.row,i.warn] = 50; quals[tmp.row,i.bad] = 0; 
}



##
## Plot the heatmap
##
if(ncol(exprs.miRNA) > 1){
  sepwidth = 0.01
  width = ncol(quals)
  height = 10
  if(width <= 5){ width = 5 }
  pdf(paste(output.dir,"/SampleQC_Heatmap.pdf",sep=""), width=width,height=height)
  par(oma=c(25,0,1,5))
  try(heatmap.2(quals, trace="n", Rowv=F, dendrogram="column", colsep=1:ncol(quals), rowsep=1:nrow(quals), sepwidth=c(sepwidth,sepwidth), col=c("red","orange","lightgreen"), cellnote=round(fracs,0), notecol="white", notecex=2, cexRow=2, cexCol=2, breaks=seq(0,100,by=33.3), key=F, keysize=0.5), silent=T)
  dev.off()
}




##
## Open PDF for diagnostic plots
##
pdf(paste(output.dir,"DiagnosticPlots.pdf",sep="/"), height=10, width=20)
#tiff(paste(output.dir,"DiagnosticPlots.tiff",sep="/"))


##
## plot distribution of clipped read lengths
##
tmp = melt(read.lengths); colnames(tmp) = c("sample","length","count")
tmp = tmp[1:max(which(tmp$count > 0)), ]
ggplot(tmp, aes(x=length, y=count, colour=sample)) +geom_line(alpha=0.75) +xlab("read length (nt)") +ylab("# reads")  +ggtitle("read-length distributions")
#ggplot(tmp, aes(x=as.factor(length), y=count)) +geom_violin()
#ggplot(tmp, aes(x=as.factor(length), y=count)) +geom_boxplot()


##
## plot distribution of # mapped reads per sample
##
tmp = log10(libSizes$all)
hist(tmp, breaks=seq(0,ceiling(max(tmp)), by=0.1), col="grey", border="white", xlab="log10(# mapped reads)", main="Library size (all mapped reads)", ylab="# samples")


##
## Plot the rRNA contamination
##
par(mfrow=c(1,2))
hist((mapping.stats$UniVec_contaminants / libSizes$all), breaks=seq(0,1,by=0.05), col="grey", border="white", xlim=c(0,1), main="UniVec contaminant signal",xlab="fraction contaminant reads",ylab="# samples")
hist((mapping.stats$rRNA / libSizes$all), breaks=seq(0,1,by=0.05), col="grey", border="white", xlim=c(0,1), main="rRNA signal",xlab="fraction rRNA reads",ylab="# samples")


##
## Calculate reads per million (RPM)
##
#libSize.use = libSizes$all
libSize.use = libSizes$miRNA
exprs.miRNA.rpm = t(10^6 * t(exprs.miRNA) / libSize.use)
exprs.tRNA.rpm = t(10^6 * t(exprs.tRNA) / libSize.use)
exprs.piRNA.rpm = t(10^6 * t(exprs.piRNA) / libSize.use)
exprs.gencode.rpm = t(10^6 * t(exprs.gencode) / libSize.use)
exprs.exogenous_miRNA.rpm = t(10^6 * t(exprs.exogenous_miRNA) / libSize.use)

par(mfrow=c(1,2), oma=c(15,0,0,0))
boxplot(log10(exprs.miRNA+1E-1), las=2, ylab="log10(miRNA counts)", main="miRNA read count")
boxplot(log10(exprs.miRNA.rpm+1E-1), las=2, ylab="log10(miRNA RPM)", main="miRNA RPM")


## Plot miRNA expression distributions
tmp = melt(log10(exprs.miRNA.rpm+1E-1))
colnames(tmp) = c("miRNA","sample","abundance")
ggplot(tmp, aes(y=abundance, x=sample, colour=sample)) +geom_violin() +geom_boxplot(alpha=0.2) +ylab("log10(RPM+1E-1)") +ggtitle("miRNA abundance distributions (RPM)") +theme(axis.ticks = element_blank(), axis.text.x = element_blank())
ggplot(tmp, aes(x=abundance, colour=sample)) +geom_density() +xlab("log10(RPM)") +ggtitle("miRNA abundance distributions (RPM)")
#ggplot(tmp, aes(y=abundance, x=sample, colour=sample)) +geom_boxplot() +ylab("log10(RPM+1E-1)") +ggtitle("miRNA abundance distributions")
#ggplot(tmp, aes(y=abundance, x=sample, colour=sample)) +geom_violin() +ylab("log10(RPM+1E-1)") +ggtitle("miRNA abundance distributions")
dev.off()


##
## Save the RPM normalised data
##
save(exprs.miRNA.rpm, exprs.tRNA.rpm, exprs.piRNA.rpm, exprs.gencode.rpm, exprs.exogenous_miRNA.rpm, file=paste(output.dir, "smallRNAQuants_ReadsPerMillion.RData", sep="/"))
write.table(exprs.miRNA.rpm, file=paste(output.dir, "miRNA_Quantifications_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)

