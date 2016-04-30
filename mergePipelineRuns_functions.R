##########################################################################################
##                                                                                      ##
## Functions to combine pipeline runs for individual samples into something more useful ##
##                                                                                      ##
## Author: Rob Kitchen (rob.kitchen@yale.edu)                                           ##
##                                                                                      ##
## Version 3.2.1 (2015-11-03)                                                           ##
##                                                                                      ##
##########################################################################################


##
## check dependencies
##
if(!"plyr" %in% rownames(installed.packages())) { install.packages("plyr",repos='http://cran.us.r-project.org') }
if(!"gplots" %in% rownames(installed.packages())) { install.packages("gplots",repos='http://cran.us.r-project.org') }
if(!"marray" %in% rownames(installed.packages())) { source("http://bioconductor.org/biocLite.R"); biocLite("marray") }
if(!"reshape2" %in% rownames(installed.packages())) { install.packages("reshape2",repos='http://cran.us.r-project.org') }
if(!"ggplot2" %in% rownames(installed.packages())) { install.packages("ggplot2",repos='http://cran.us.r-project.org') }
if(!"tools" %in% rownames(installed.packages())) { install.packages("tools",repos='http://cran.us.r-project.org') }
require(plyr)
require(gplots)
require(marray)
require(reshape2)
require(ggplot2)
require(tools)



##
## Function to recursively search a given directory for pipeline output
##
SearchForSampleData = function(base.dir, directory=""){
  to.return = NULL
  dir.use = paste(base.dir, directory, sep = "/")
  subdirs = dir(dir.use)
  if(length(subdirs) > 0){
    i.stats = grep("\\.stats$", subdirs, perl=T)
    i.zip = grep("\\.zip$", subdirs, perl=T)
    i.tar = grep("\\.tgz$|\\.tar.gz$", subdirs, perl=T)
    
    ## handle decompressed pipeline output
    if(length(i.stats) > 0){
      tmp = gsub("\\.stats$","",subdirs[i.stats])
      to.return = c(to.return, paste(dir.use,tmp,sep="/"))
    }
    
    ## handle zipped pipeline output
    if(length(i.zip) > 0){
      for(x in i.zip){
        tmp.contents = unzip(paste(dir.use,subdirs[x],sep="/"), list=T)[,1]
        if(length(grep("\\.stats$", tmp.contents, perl=T)) > 0){
          try(unzip(paste(dir.use,subdirs[x],sep="/"), exdir=gsub("\\.zip","",file_path_as_absolute(paste(dir.use,subdirs[x],sep="/"))), overwrite=FALSE), silent=T)
          to.return = c(to.return, paste(dir.use,gsub("\\.zip$","",subdirs[x]),sep="/",gsub("\\.stats$","",tmp.contents[grep("\\.stats$", tmp.contents)])))
        }
      }
    }
    
    ## handle [tar] gzipped pipeline output
    if(length(i.tar) > 0){
      for(x in i.tar){
        tmp.dir = paste(dir.use,subdirs[x],sep="/")
        tmp.contents = untar(tmp.dir, list=T, tar="tar")
        if(length(grep("\\.stats$", tmp.contents, perl=T)) > 0){
          try(untar(tmp.dir, exdir=gsub("\\.tgz$|\\.tar.gz$","",file_path_as_absolute(tmp.dir)), tar="tar"), silent=T)
          to.return = c(to.return, paste(dir.use, gsub("\\.tgz$|\\.tar.gz$","",subdirs[x]), gsub("\\.stats$","",tmp.contents[grep("\\.stats$", tmp.contents)]),sep="/"))
        }
      }
    }
    
    ## handle unknown directories
    i.known = c(i.stats,i.zip,i.tar)
    if(length(i.known) == 0){		# there are no .stats, .zip, or .tgz/.tar.gz files in the directory!
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
## Prints the given message with a timestamp
##
printMessage = function(message=""){
  cat(as.character(Sys.time()),":  ",paste(message,sep=""),"\n",sep="")
}



##
##
##
processSamplesInDir = function(data.dir, output.dir=data.dir){
  
  ##  Look for samples to merge
  printMessage("Searching for valid exceRpt pipeline output...")
  samplePathList = unique(SearchForSampleData(data.dir,""))
  
  ## get sample names and remove duplicates:
  sampleIDs = sapply(samplePathList, function(path){ tmp=unlist(strsplit(path,"/")); tmp[length(tmp)] })
  samplePathList = samplePathList[!duplicated(sampleIDs)]
  
  ## -- Kill script if we do not have any samples to process
  NumberOfCompatibleSamples = length(samplePathList)
  stopifnot(NumberOfCompatibleSamples > 0)
  printMessage(c("Found ",NumberOfCompatibleSamples," valid samples"))
  
  
  ##
  ## Create objects to contain the data
  ##
  sample.data = vector(mode="list",length=length(samplePathList))
  allIDs.miRNA = NULL
  allIDs.tRNA = NULL
  allIDs.piRNA = NULL
  allIDs.gencode = NULL
  allIDs.circularRNA = NULL
  allIDs.exogenous_miRNA = NULL
  allIDs.exogenous_genomes = NULL
  mapping.stats = matrix(0,nrow=length(samplePathList),ncol=30, dimnames=list(1:length(samplePathList), c("input","successfully_clipped","failed_quality_filter","failed_homopolymer_filter","calibrator","UniVec_contaminants","rRNA","reads_used_for_alignment","genome","miRNA_sense","miRNA_antisense","miRNAprecursor_sense","miRNAprecursor_antisense","tRNA_sense","tRNA_antisense","piRNA_sense","piRNA_antisense","gencode_sense","gencode_antisense","circularRNA_sense","circularRNA_antisense","not_mapped_to_genome_or_libs","repetitiveElements","endogenous_gapped","input_to_exogenous_miRNA","exogenous_miRNA","input_to_exogenous_rRNA","exogenous_rRNA","input_to_exogenous_genomes","exogenous_genomes")))
  maxReadLength = 1000
  read.lengths = matrix(0,nrow=length(samplePathList),ncol=maxReadLength+1,dimnames=list(1:length(samplePathList), 0:maxReadLength))
  
  
  ##
  ## Loop through all samples and read the pipeline output
  ##
  printMessage(c("Reading sample data..."))
  removeSamples = NULL
  for(i in 1:length(samplePathList)){
    ## Parse the sampleID from the path:
    tmp = unlist(strsplit(samplePathList[i], "/"))
    thisSampleID = tmp[length(tmp)]
    
    ## Get timings and check this sample finished successfully
    tmp.stats = read.table(paste(samplePathList[i],".stats",sep=""), stringsAsFactors=F, fill=T, header=F, sep="\t",skip=0,comment.char="")
    x.start = grep("#STATS",tmp.stats[,1])
    x.end = grep("#END OF STATS",tmp.stats[,1])
    if(length(x.start) > 0  &&  length(x.end) > 0){
      tmp.start = strptime(unlist(strsplit(tmp.stats[x.start[1],1],"Run started at "))[2],"%Y-%m-%d--%H:%M:%S")
      tmp.end = strptime(unlist(strsplit(tmp.stats[x.end[1],1],"Run completed at "))[2],"%Y-%m-%d--%H:%M:%S")
      runTiming = data.frame(start=tmp.start, completed=tmp.end, duration=difftime(tmp.end,tmp.start), duration_secs=as.numeric(difftime(tmp.end,tmp.start,units="secs")))
      continue = T
    }else{
      continue = F
      removeSamples = c(removeSamples, i)
      printMessage(c("[",i,"/",length(samplePathList),"] WARNING: Incomplete run for sample \'",thisSampleID,"\', ignoring"))
    }
    
    if(continue == T){
      ## Read sample mapping stats
      tmp.stats = read.table(paste(samplePathList[i],".stats",sep=""), stringsAsFactors=F, fill=T, header=T, sep="\t",skip=0)
      tmp.stats[tmp.stats[,1] %in% "clipped", 1] = "successfully_clipped"
      #mapping.stats[i, match(tmp.stats[,1], colnames(mapping.stats))] = as.numeric(tmp.stats[,2])
      mapping.stats[i, match(tmp.stats[,1], colnames(mapping.stats))] = as.numeric(tmp.stats[,2])
      rownames(mapping.stats)[i] = thisSampleID
      
      ##
      ## Read the clipped read lengths  
      ##
      if(length(grep(".readLengths.txt$", dir(samplePathList[i]))) == 1){
        tmp = read.table(paste(samplePathList[i], dir(samplePathList[i])[grep(".readLengths.txt$", dir(samplePathList[i]))], sep="/"))
        read.lengths[i, 1:ncol(tmp)] = as.numeric(tmp[2,])
        rownames(read.lengths)[i] = thisSampleID
      }
      
      
      ##
      ## Read the adapter sequence
      ##
      if(paste(thisSampleID,".adapterSeq",sep="") %in% dir(samplePathList[i])){
        tmp.seq = try(read.table(paste(samplePathList[i],"/",thisSampleID,".adapterSeq",sep="")), silent=T)
        if(class(tmp.seq) == "try-error"){
          adapterSeq = NA
        }else{
          adapterSeq = as.character(tmp.seq[1,1])
        }
      }
      
      
      ##
      ## Read sample data
      ##
      availableFiles = dir(samplePathList[i])
      miRNA_sense=miRNA_antisense = tRNA_sense=tRNA_antisense = piRNA_sense=piRNA_antisense = gencode_sense=gencode_antisense = circRNA_sense=circRNA_antisense = NULL
      
      if("readCounts_miRNAmature_sense.txt" %in% availableFiles){
        miRNA_sense = read.table(paste(samplePathList[i],"readCounts_miRNAmature_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        miRNA_sense = cbind(miRNA_sense, ID=sapply(rownames(miRNA_sense), function(id){ multiID = unlist(strsplit(id,"\\|")); multiIDs = sapply(multiID, function(idPart){unlist(strsplit(idPart,":"))[1]});   if(length(multiIDs) == 1){ multiIDs }else{ paste(sort(multiIDs),collapse="|") }  }))
      }
      if("readCounts_miRNAmature_antisense.txt" %in% availableFiles){
        miRNA_antisense = read.table(paste(samplePathList[i],"readCounts_miRNAmature_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        miRNA_antisense = cbind(miRNA_antisense, ID=sapply(rownames(miRNA_antisense), function(id){ multiID = unlist(strsplit(id,"\\|")); multiIDs = sapply(multiID, function(idPart){unlist(strsplit(idPart,":"))[1]});   if(length(multiIDs) == 1){ multiIDs }else{ paste(sort(multiIDs),collapse="|") }  }))
      }
      
      if("readCounts_tRNA_sense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_tRNA_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), function(id){ unlist(strsplit(id,"-"))[2]  }))
        tRNA_sense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(tRNA_sense)[-1] = colnames(tmp)[1:4]
        tRNA_sense = tRNA_sense[order(tRNA_sense$multimapAdjustedReadCount,decreasing=T), ]
      }
      if("readCounts_tRNA_antisense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_tRNA_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), function(id){ unlist(strsplit(id,"-"))[2]  }))
        tRNA_antisense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(tRNA_antisense)[-1] = colnames(tmp)[1:4]
        tRNA_antisense = tRNA_antisense[order(tRNA_antisense$multimapAdjustedReadCount,decreasing=T), ]
      }
      
      if("readCounts_piRNA_sense.txt" %in% availableFiles){
        piRNA_sense = read.table(paste(samplePathList[i],"readCounts_piRNA_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      if("readCounts_piRNA_antisense.txt" %in% availableFiles){
        piRNA_antisense = read.table(paste(samplePathList[i],"readCounts_piRNA_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      
      
      makeGeneID = function(id){ 
        bits=unlist(strsplit(id,":")); geneNameBits=unlist(strsplit(bits[3],"-")); 
        geneName = geneNameBits[1]
        if(length(geneNameBits) > 2){ geneName=paste(geneNameBits[-length(geneNameBits)],collapse="-") }
        paste(geneName,bits[2],sep=":") 
      }
      if("readCounts_gencode_sense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_gencode_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), makeGeneID))
        gencode_sense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(gencode_sense)[-1] = colnames(tmp)[1:4]
        gencode_sense = gencode_sense[order(gencode_sense$multimapAdjustedReadCount,decreasing=T), ]
      }
      if("readCounts_gencode_antisense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_gencode_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), makeGeneID))
        gencode_antisense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(gencode_antisense)[-1] = colnames(tmp)[1:4]
        gencode_antisense = gencode_antisense[order(gencode_antisense$multimapAdjustedReadCount,decreasing=T), ]
      }
      
      if("readCounts_circRNA_sense.txt" %in% availableFiles){
        circRNA_sense = read.table(paste(samplePathList[i],"readCounts_circRNA_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      if("readCounts_circRNA_antisense.txt" %in% availableFiles){
        circRNA_antisense = read.table(paste(samplePathList[i],"readCounts_circRNA_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      
      
      # Update list of detected smallRNA IDs
      allIDs.miRNA = unique(c(allIDs.miRNA, as.character(miRNA_sense$ID)))
      allIDs.tRNA = unique(c(allIDs.tRNA, as.character(tRNA_sense$ID)))
      allIDs.piRNA = unique(c(allIDs.piRNA, rownames(piRNA_sense)))
      allIDs.gencode = unique(c(allIDs.gencode, as.character(gencode_sense$ID)))
      allIDs.circularRNA = unique(c(allIDs.circularRNA, rownames(circRNA_sense)))
      
      sample.data[[i]] = list("miRNA_sense"=miRNA_sense,"miRNA_antisense"=miRNA_antisense, "tRNA_sense"=tRNA_sense,"tRNA_antisense"=tRNA_antisense, "piRNA_sense"=piRNA_sense,"piRNA_antisense"=piRNA_antisense, "gencode_sense"=gencode_sense,"gencode_antisense"=gencode_antisense, "circRNA_sense"=circRNA_sense,"circRNA_antisense"=circRNA_antisense, "adapterSeq"=adapterSeq, "runTiming"=runTiming)
      names(sample.data)[i] = thisSampleID
      
      
      ##
      ## Read exogenous miRNA alignments (if applicable)
      ##
      if("EXOGENOUS_miRNA" %in% availableFiles){
        tmp.dir = paste(samplePathList[i],"EXOGENOUS_miRNA",sep="/")
        if("mature_sense.grouped" %in% dir(tmp.dir)){
          #sample.data[[i]]$exogenous_miRNA = readData(tmp.dir,"mature_sense.grouped", 4)
          #allIDs.exogenous_miRNA = unique(c(allIDs.exogenous_miRNA, as.character(sample.data[[i]]$exogenous_miRNA[,1])))
        }
      }
      
      
      ##
      ## Read exogenous genome alignments (if applicable)
      ##
      if("EXOGENOUS_genomes" %in% availableFiles){
        tmp.dir = paste(samplePathList[i],"EXOGENOUS_genomes",sep="/")
        if("ExogenousGenomicAlignments.result.txt" %in% dir(tmp.dir)){
          tmp = read.table(paste(tmp.dir, "ExogenousGenomicAlignments.result.txt",sep="/"), stringsAsFactors=F, comment.char="",header=T)
          rownames(tmp) = apply(tmp[,1:2], 1, paste, collapse="|")
          #colnames(tmp) = c("Kingdom","Species","ReadCount_all","ReadCount_kingdomSpecific","ReadCount_speciesSpecific")
          sample.data[[i]]$exogenous_genomes = tmp
          allIDs.exogenous_genomes = unique(c(allIDs.exogenous_genomes, as.character(rownames(sample.data[[i]]$exogenous_genomes))))
        }
      }
      printMessage(c("[",i,"/",length(samplePathList),"] Added sample \'",thisSampleID,"\'"))
    }
  }
  
  
  
  ##
  ## Remove failed/incomplete samples
  ##
  stopifnot(length(removeSamples) < length(sample.data))
  if(length(removeSamples) > 0){
    read.lengths = read.lengths[-removeSamples, ]
    sample.data = sample.data[-removeSamples]
    mapping.stats = mapping.stats[-removeSamples,]
  }
  
  
  
  ##
  ## Trim read-length matrix
  ##
  read.lengths = read.lengths[,0:(max(as.numeric(colnames(read.lengths[, colSums(read.lengths) > 0])))+1)]
  #read.lengths = read.lengths[,colSums(read.lengths) > 0, drop=F]
  
  
  ##
  ## Collect IDs
  ##
  allIDs = list("miRNA_sense"=allIDs.miRNA, "tRNA_sense"=allIDs.tRNA, "piRNA_sense"=allIDs.piRNA, "gencode_sense"=allIDs.gencode, "circRNA_sense"=allIDs.circularRNA, "exogenous_miRNA"=allIDs.exogenous_miRNA, "exogenous_genomes"=allIDs.exogenous_genomes)
  
  
  ##
  ## Convert sample data to large per-smallRNA expression matrices
  ##
  printMessage("Creating raw read-count matrices for available libraries")
  #run.duration = data.frame(runDuration_string=rep("",length(sample.data)), runDuration_secs=rep(0,length(sample.data)),stringsAsFactors = F)
  run.duration = data.frame(runDuration_secs=rep(0,length(sample.data)),stringsAsFactors = F)
  rownames(run.duration) = names(sample.data)
  exprs.miRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$miRNA_sense), dimnames=list(allIDs$miRNA_sense, names(sample.data)))
  exprs.tRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$tRNA_sense), dimnames=list(allIDs$tRNA_sense, names(sample.data)))
  exprs.piRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$piRNA_sense), dimnames=list(allIDs$piRNA_sense, names(sample.data)))
  exprs.gencode = matrix(0,ncol=length(sample.data),nrow=length(allIDs$gencode_sense), dimnames=list(allIDs$gencode_sense, names(sample.data)))
  exprs.circRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$circRNA_sense), dimnames=list(allIDs$circRNA_sense, names(sample.data)))
  #exprs.exogenous_miRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$exogenous_miRNA), dimnames=list(allIDs$exogenous_miRNA, names(sample.data)))
  exprs.exogenousGenomes_speciesSpecific = matrix(0,ncol=length(sample.data),nrow=length(allIDs$exogenous_genomes), dimnames=list(allIDs$exogenous_genomes, names(sample.data)))
  exprs.exogenousGenomes_kingdomSpecific = matrix(0,ncol=length(sample.data),nrow=length(allIDs$exogenous_genomes), dimnames=list(allIDs$exogenous_genomes, names(sample.data)))
  for(i in 1:length(sample.data)){
    run.duration[i,] = sample.data[[i]]$runTiming[1,4,drop=F]
    exprs.miRNA[match(sample.data[[i]]$miRNA_sense$ID, rownames(exprs.miRNA)),i] = as.numeric(sample.data[[i]]$miRNA_sense$multimapAdjustedReadCount)
    exprs.tRNA[match(sample.data[[i]]$tRNA_sense$ID, rownames(exprs.tRNA)),i] = as.numeric(sample.data[[i]]$tRNA_sense$multimapAdjustedReadCount)
    exprs.piRNA[match(rownames(sample.data[[i]]$piRNA_sense), rownames(exprs.piRNA)),i] = as.numeric(sample.data[[i]]$piRNA_sense$multimapAdjustedReadCount)
    exprs.gencode[match(sample.data[[i]]$gencode_sense$ID, rownames(exprs.gencode)),i] = as.numeric(sample.data[[i]]$gencode_sense$multimapAdjustedReadCount)
    exprs.circRNA[match(rownames(sample.data[[i]]$circRNA_sense), rownames(exprs.circRNA)),i] = as.numeric(sample.data[[i]]$circRNA_sense$multimapAdjustedReadCount)
    #exprs.exogenous_miRNA[match(sample.data[[i]]$exogenous_miRNA[,1], rownames(exprs.exogenous_miRNA)),i] = as.numeric(sample.data[[i]]$exogenous_miRNA[,2])
    exprs.exogenousGenomes_speciesSpecific[match(rownames(sample.data[[i]]$exogenous_genomes), rownames(exprs.exogenousGenomes_speciesSpecific)),i] = as.numeric(sample.data[[i]]$exogenous_genomes[,5])
    exprs.exogenousGenomes_kingdomSpecific[match(rownames(sample.data[[i]]$exogenous_genomes), rownames(exprs.exogenousGenomes_kingdomSpecific)),i] = as.numeric(sample.data[[i]]$exogenous_genomes[,4])
  }
  
  
  ##
  ## Calculate the total number of mapped reads to the rRNA, genome, and exogenous sequences
  ##
  mapping.stats[is.na(mapping.stats)] = 0
  mapping.stats = as.data.frame(mapping.stats)
  libSizes = list()
  libSizes$input = mapping.stats[,colnames(mapping.stats) %in% c("input")]
  libSizes$all = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("rRNA","genome","miRNA_exogenous_sense")])
  libSizes$endogenous = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("rRNA","genome")])
  libSizes$genome = mapping.stats[,colnames(mapping.stats) %in% "genome"]
  libSizes$smRNA = mapping.stats[,grep("sense",colnames(mapping.stats))]
  libSizes$miRNA = colSums(exprs.miRNA)
  
  
  ##
  ## Save the raw count data
  ##
  printMessage("Saving raw data to disk")
  #save(exprs.miRNA, exprs.tRNA, exprs.piRNA, exprs.gencode, exprs.circRNA, exprs.exogenous_miRNA, exprs.exogenous_genomes, mapping.stats, libSizes, read.lengths, file=paste(output.dir, "exceRpt_smallRNAQuants_ReadCounts.RData", sep="/"))
  save(exprs.miRNA, exprs.tRNA, exprs.piRNA, exprs.gencode, exprs.circRNA, exprs.exogenousGenomes_speciesSpecific, exprs.exogenousGenomes_kingdomSpecific, mapping.stats, libSizes, read.lengths, file=paste(output.dir, "exceRpt_smallRNAQuants_ReadCounts.RData", sep="/"))
  write.table(exprs.miRNA, file=paste(output.dir, "exceRpt_miRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.tRNA, file=paste(output.dir, "exceRpt_tRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.piRNA, file=paste(output.dir, "exceRpt_piRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.gencode, file=paste(output.dir, "exceRpt_gencode_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.circRNA, file=paste(output.dir, "exceRpt_circularRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.exogenousGenomes_speciesSpecific, file=paste(output.dir, "exceRpt_exogenousGenomes_speciesSpecific_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.exogenousGenomes_kingdomSpecific, file=paste(output.dir, "exceRpt_exogenousGenomes_kingdomSpecific_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(read.lengths, file=paste(output.dir, "exceRpt_ReadLengths.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  
  ##
  ## Calculate the fractions of reads mapping at each stage
  ##
  write.table(mapping.stats, file=paste(output.dir,"exceRpt_readMappingSummary.txt",sep="/"), sep="\t", col.names=NA, quote=F)
  
  
  ##
  ## Order samples based on similarity of mapping statistics
  ##
  sampleOrder = 1
  if(nrow(mapping.stats) > 1){
    h = hclust(dist(1-cor(t(mapping.stats))))
    sampleOrder = h$order
  }
  
  
  ##
  ## Guess sample grouping based on their names
  ##
  #sampleNames = rownames(mapping.stats)
  #tmp = adist(sampleNames,ignore.case=T)
  #dimnames(tmp) = list(sampleNames,sampleNames)
  #par(mfrow=c(1,1))
  #plot(hclust(as.dist(tmp)))
  
  ##
  ## Calculate optimal number of sample groups
  ##
  #require(cluster)
  #K.max=10
  #B=200
  ## Slightly faster way to use pam (see below)
  #pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
  #gs.kmeans <- clusGap(t(tmp), FUN=kmeans, nstart=3, K.max=K.max, B=B)
  #gs.pam <- clusGap(t(tmp), FUN=pam1, K.max=K.max, B=B)
  ## plot?
  ##par(mfrow=c(2,1))
  ##ylim=c(0.3,0.45)
  ##ylim=NULL
  ##plot(gs.kmeans, ylim=ylim, main="kmeans"); abline(h=0.398, col="darkgreen")
  ##plot(gs.pam, ylim=ylim, main="pam"); abline(h=0.406, col="darkgreen")
  ## find optimal cluster number
  #nGroups.kmeans = as.numeric(names(sort(table(sapply(c(1/4, 1,2,4), function(SEf){ sapply(eval(formals(maxSE)$method), function(M){maxSE(gs.kmeans$Tab[,3], gs.kmeans$Tab[,4], method = M, SE.factor = SEf)}) })),decreasing=T)[1]))
  #nGroups.pam = as.numeric(names(sort(table(sapply(c(1/4, 1,2,4), function(SEf){ sapply(eval(formals(maxSE)$method), function(M){maxSE(gs.pam$Tab[,3], gs.pam$Tab[,4], method=M, SE.factor=SEf)}) })),decreasing=T)[1]))
  #as.numeric(names(sort(table(sapply(c(1/4, 1,2,4), function(SEf){ sapply(eval(formals(maxSE)$method), function(M){maxSE(gs.kmeans$Tab[,3], gs.kmeans$Tab[,4], method = M, SE.factor = SEf)}) })[3,]),decreasing=T)[1]))
  #as.numeric(names(sort(table(sapply(c(1/4, 1,2,4), function(SEf){ sapply(eval(formals(maxSE)$method), function(M){maxSE(gs.pam$Tab[,3], gs.pam$Tab[,4], method = M, SE.factor = SEf)}) })[3,]),decreasing=T)[1]))
  ##if(nGroups.pam <= 5){
  ## 
  #  
  #}
  
  
  ##
  ## Open PDF for diagnostic plots
  ##
  printMessage("Creating QC plots")
  pdf(paste(output.dir,"exceRpt_DiagnosticPlots.pdf",sep="/"), height=10, width=20)
  #tiff(paste(output.dir,"DiagnosticPlots.tiff",sep="/"))
  
  
  if(ncol(read.lengths) > 1){
    ##
    ## plot distribution of clipped read lengths - read count
    ##
    tmp = melt(read.lengths); colnames(tmp) = c("sample","length","count")
    #tmp = tmp[1:max(which(tmp$count > 0)), ]
    p = ggplot(tmp, aes(x=length, y=count, colour=sample)) +geom_line(alpha=0.75) +xlab("read length (nt)") +ylab("# reads") +ggtitle("read-length distributions") +xlim(14,min(c(75,max(tmp$length))))
    if(nrow(read.lengths) > 30){ p = p +guides(colour=FALSE) }
    print(p)
    #ggplot(tmp, aes(x=as.factor(length), y=count)) +geom_violin()
    #ggplot(tmp, aes(x=as.factor(length), y=count)) +geom_boxplot()
    
    
    ##
    ## plot distribution of clipped read lengths - fraction
    ##
    tmp = melt(t(apply(read.lengths, 1, function(row){ row/sum(row) }))); colnames(tmp) = c("sample","length","fraction")
    #tmp = tmp[1:max(which(tmp$fraction > 0)), ]
    p = ggplot(tmp, aes(x=length, y=fraction, colour=sample)) +geom_line(alpha=0.75) +xlab("read length (nt)") +ylab("fraction of reads") +ggtitle("read-length distributions") +xlim(14,min(c(75,max(tmp$length))))
    if(nrow(read.lengths) > 30){ p = p +guides(colour=FALSE) }
    print(p)
  }
  
  
  
  ##
  ## Plot run duration of each sample
  ##
  tmp=melt(as.matrix(run.duration))
  colnames(tmp) = c("sampleID","stuff","runDuration_seconds")
  tmp = cbind(tmp, category=.bincode(tmp[,3], breaks=c(0,as.numeric(quantile(tmp[,3],probs=c(0.10,0.90,1))))))
  tmp = cbind(tmp, colour=tmp$category)
  tmp = cbind(tmp, inputReadCount=mapping.stats$input)
  tmp$category[tmp$category == 1] = "fast"
  tmp$category[tmp$category == 2] = "normal"
  tmp$category[tmp$category == 3] = "slow"
  tmp$colour[tmp$colour == 1] = "red"
  tmp$colour[tmp$colour == 2] = "green"
  tmp$colour[tmp$colour == 3] = "blue"
  p = ggplot(tmp, aes(x=sampleID,y=runDuration_seconds,fill=colour)) +geom_bar(stat="identity") +facet_grid(~category,scales="free_x",space="free_x") +guides(fill=FALSE) +theme(axis.text.x=element_text(angle=90, hjust=1,vjust=0.5))
  print(p)
  
  p = ggplot(tmp, aes(x=inputReadCount,y=runDuration_seconds,colour=colour)) +geom_point(size=10) +guides(colour=FALSE) +scale_y_log10(limits=c(1,10^ceiling(log10(max(tmp$runDuration_seconds)))), breaks=10^seq(1:ceiling(log10(max(tmp$runDuration_seconds))))) +scale_x_log10(limits=c(min(c(100000,10^floor(log10(min(tmp$inputReadCount))))),10^ceiling(log10(max(tmp$inputReadCount)))), breaks=10^seq(min(c(100000,floor(log10(min(tmp$inputReadCount))))),ceiling(log10(max(tmp$inputReadCount)))))
  print(p)
  
  
  ##
  ## plot distribution of # mapped reads per sample
  ##
  tmp = log10(libSizes$all)
  hist(tmp, breaks=seq(0,ceiling(max(tmp)), by=0.1), col="grey", border="white", xlab="log10(# mapped reads)", main="Library size (all mapped reads)", ylab="# samples")
  
  
  ##
  ## Plot the rRNA contamination
  ##
  #par(mfrow=c(1,2))
  #hist((mapping.stats$UniVec_contaminants / libSizes$all), breaks=seq(0,1,by=0.05), col="grey", border="white", xlim=c(0,1), main="UniVec contaminant signal",xlab="fraction contaminant reads",ylab="# samples")
  #hist((mapping.stats$rRNA / libSizes$all), breaks=seq(0,1,by=0.05), col="grey", border="white", xlim=c(0,1), main="rRNA signal",xlab="fraction rRNA reads",ylab="# samples")
  
  
  
  mapping.stats.orig = mapping.stats
  mapping.stats = mapping.stats[,-grep("input_to_",colnames(mapping.stats))]
  
  
  
  ##
  ## Plot heatmap of mapping percentages through the pipeline
  ##
  toplot = melt(as.matrix(mapping.stats / mapping.stats$input)); colnames(toplot) = c("Sample","Stage","ReadFraction")
  toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
  toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
  p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +ggtitle("fraction aligned reads (normalised by # input reads)")
  if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
  print(p)
  
  ##
  ## Plot heatmap of mapping percentages through the pipeline
  ##
  if(max(mapping.stats$successfully_clipped) > 0){
    toplot = melt(as.matrix(mapping.stats / mapping.stats$successfully_clipped)[,-1,drop=F]); colnames(toplot) = c("Sample","Stage","ReadFraction")
    toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
    toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
    p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +ggtitle("fraction aligned reads (normalised by # adapter-clipped reads)")
    if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
    print(p)
  }
  
  ##
  ## Plot heatmap of mapping percentages through the pipeline
  ##
  toplot = melt(as.matrix(mapping.stats / mapping.stats$reads_used_for_alignment)[,-c(1:7),drop=F]); colnames(toplot) = c("Sample","Stage","ReadFraction")
  toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
  toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
  p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +ggtitle("fraction aligned reads (normalised by # non-contaminant reads)")
  if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
  print(p)
  
  
  
  ##
  ## Plot breakdown of counts in each biotype 
  ##
  require(plyr)
  sampleTotals=matrix(NA,ncol=nrow(mapping.stats),nrow=0); colnames(sampleTotals) = rownames(mapping.stats)
  if(nrow(exprs.miRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.miRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "miRNA"
  }
  if(nrow(exprs.tRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.tRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "tRNA"
  }
  if(nrow(exprs.piRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.piRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "piRNA"
  }
  if(nrow(exprs.gencode) > 0){
    tmp = data.frame(biotype=sapply(rownames(exprs.gencode), function(id){bits=unlist(strsplit(id,":")); bits[length(bits)]}), exprs.gencode)
    tmp = ddply(tmp, "biotype", function(mat){ colSums(mat[,-1,drop=F]) })
    rownames(tmp) = tmp[,1]; tmp = tmp[,-1,drop=F]
    colnames(tmp) = colnames(sampleTotals)
    ## add gencode miRNA to the existing count...
    if("miRNA" %in% rownames(tmp)  &&  "miRNA" %in% rownames(sampleTotals)){
      i = which(rownames(tmp) %in% "miRNA")
      j = which(rownames(sampleTotals) %in% "miRNA")
      for(x in 1:ncol(sampleTotals)){
        sampleTotals[j,x] = sampleTotals[j,x] + tmp[i,x]
      }
      tmp = tmp[-i,,drop=F]
    }
    sampleTotals = rbind(sampleTotals, tmp)
  }
  if(nrow(exprs.circRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.circRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "circularRNA"
  }
  sampleTotals = rbind(sampleTotals, mapping.stats$exogenous_miRNA)
  rownames(sampleTotals)[nrow(sampleTotals)] = "exogenous_miRNA"
  sampleTotals = rbind(sampleTotals, mapping.stats$exogenous_rRNA)
  rownames(sampleTotals)[nrow(sampleTotals)] = "exogenous_rRNA"
  sampleTotals = rbind(sampleTotals, mapping.stats$exogenous_genomes)
  rownames(sampleTotals)[nrow(sampleTotals)] = "exogenous_genomes"
  
  sampleTotals = sampleTotals[order(apply(sampleTotals, 1, median, na.rm=T), decreasing=F), ,drop=F]
  tmp = melt(as.matrix(sampleTotals))
  colnames(tmp) = c("biotype","sampleID","readCount")
  p = ggplot(na.omit(tmp), aes(y=readCount,x=biotype, colour=biotype)) +geom_hline(y=1,linetype="dashed") +geom_boxplot() +scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000,100000,1000000,10000000,100000000)) +guides(colour=FALSE) +coord_flip()
  print(p)
  
  
  
  ##
  ## Calculate reads per million (RPM)
  ##
  #libSize.use = libSizes$all
  #libSize.use = libSizes$miRNA
  libSize.use = libSizes$genome
  exprs.miRNA.rpm = t(10^6 * t(exprs.miRNA) / libSize.use)
  exprs.tRNA.rpm = t(10^6 * t(exprs.tRNA) / libSize.use)
  exprs.piRNA.rpm = t(10^6 * t(exprs.piRNA) / libSize.use)
  exprs.gencode.rpm = t(10^6 * t(exprs.gencode) / libSize.use)
  #exprs.repElements.rpm = t(10^6 * t(exprs.repElements) / libSize.use)
  #exprs.exogenous_miRNA.rpm = t(10^6 * t(exprs.exogenous_miRNA) / libSize.use)
  exprs.exogenousGenomes_speciesSpecific.rpm = t(10^6 * t(exprs.exogenousGenomes_speciesSpecific) / libSize.use)
  exprs.exogenousGenomes_kingdomSpecific.rpm = t(10^6 * t(exprs.exogenousGenomes_kingdomSpecific) / libSize.use)
  
  #par(mfrow=c(1,2), oma=c(15,0,0,0))
  #boxplot(log10(exprs.miRNA[, sampleOrder]+1E-1), las=2, ylab="log10(miRNA counts)", main="miRNA read count")
  #boxplot(log10(exprs.miRNA.rpm[, sampleOrder]+1E-1), las=2, ylab="log10(miRNA RPM)", main="miRNA RPM")
  
  
  ## Plot miRNA expression distributions
  if(nrow(exprs.miRNA) > 0){
    tmp = melt(exprs.miRNA)
    colnames(tmp) = c("miRNA","sample","abundance")
    p = ggplot(tmp, aes(y=abundance, x=sample, colour=sample)) +geom_violin() +geom_boxplot(alpha=0.2) +ylab("Read count") +ggtitle("miRNA abundance distributions (raw counts)") +theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +scale_y_log10()
    if(ncol(exprs.miRNA.rpm) > 30){ p = p +guides(colour=FALSE) }
    print(p)
    
    p = ggplot(tmp, aes(x=abundance, colour=sample)) +geom_density() +xlab("Read count") +ggtitle("miRNA abundance distributions (raw counts)") +scale_x_log10()
    if(ncol(exprs.miRNA.rpm) > 30){ p = p +guides(colour=FALSE) }
    print(p)
    
    tmp = melt(exprs.miRNA.rpm)
    colnames(tmp) = c("miRNA","sample","abundance")
    p = ggplot(tmp, aes(y=abundance, x=sample, colour=sample)) +geom_violin() +geom_boxplot(alpha=0.2) +ylab("Reads per million (RPM)") +ggtitle("miRNA abundance distributions (RPM)") +theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +scale_y_log10()
    if(ncol(exprs.miRNA.rpm) > 30){ p = p +guides(colour=FALSE) }
    print(p)
    
    p = ggplot(tmp, aes(x=abundance, colour=sample)) +geom_density() +xlab("Reads per million (RPM)") +ggtitle("miRNA abundance distributions (RPM)") +scale_x_log10()
    if(ncol(exprs.miRNA.rpm) > 30){ p = p +guides(colour=FALSE) }
    print(p)
  }
  
  
  
  ##
  ## Finally, plot exogenous if there are any
  ##
  if(nrow(exprs.exogenousGenomes_speciesSpecific) > 0  &&  ncol(exprs.exogenousGenomes_speciesSpecific) > 1){
    par(oma=c(5,0,0,8))
    tmp.order = order(apply(t(t(exprs.exogenousGenomes_speciesSpecific)/colSums(exprs.exogenousGenomes_speciesSpecific)), 1, median), decreasing=T)
    heatmap.2(log10(exprs.exogenousGenomes_speciesSpecific[tmp.order, ][1:100,]+0.1),trace="none")
    tmp.order = order(apply(t(t(exprs.exogenousGenomes_kingdomSpecific)/colSums(exprs.exogenousGenomes_kingdomSpecific)), 1, median), decreasing=T)
    heatmap.2(log10(exprs.exogenousGenomes_kingdomSpecific[tmp.order, ][1:100,]+0.1),trace="none")
    
    tmp.order = order(apply(t(t(exprs.exogenousGenomes_speciesSpecific.rpm)/colSums(exprs.exogenousGenomes_speciesSpecific.rpm)), 1, median), decreasing=T)
    heatmap.2(log10(exprs.exogenousGenomes_speciesSpecific.rpm[tmp.order, ][1:100,]+0.1),trace="none")
    tmp.order = order(apply(t(t(exprs.exogenousGenomes_kingdomSpecific.rpm)/colSums(exprs.exogenousGenomes_kingdomSpecific.rpm)), 1, median), decreasing=T)
    heatmap.2(log10(exprs.exogenousGenomes_kingdomSpecific.rpm[tmp.order, ][1:100,]+0.1),trace="none") 
  }
  dev.off()
  
  
  ##
  ## Save the RPM normalised data
  ##
  printMessage("Saving normalised data to disk")
  #save(exprs.miRNA.rpm, exprs.tRNA.rpm, exprs.piRNA.rpm, exprs.gencode.rpm, exprs.exogenous_miRNA.rpm, exprs.exogenous_genomes.rpm, file=paste(output.dir, "exceRpt_smallRNAQuants_ReadsPerMillion.RData", sep="/"))
  save(exprs.miRNA.rpm, exprs.tRNA.rpm, exprs.piRNA.rpm, exprs.gencode.rpm, exprs.exogenousGenomes_speciesSpecific.rpm, exprs.exogenousGenomes_kingdomSpecific.rpm, file=paste(output.dir, "exceRpt_smallRNAQuants_ReadsPerMillion.RData", sep="/"))
  write.table(exprs.miRNA.rpm, file=paste(output.dir, "exceRpt_miRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.tRNA.rpm, file=paste(output.dir, "exceRpt_tRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.piRNA.rpm, file=paste(output.dir, "exceRpt_piRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.gencode.rpm, file=paste(output.dir, "exceRpt_gencode_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.exogenousGenomes_speciesSpecific.rpm, file=paste(output.dir, "exceRpt_exogenousGenomes_speciesSpecific_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(exprs.exogenousGenomes_kingdomSpecific.rpm, file=paste(output.dir, "exceRpt_exogenousGenomes_kingdomSpecific_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  printMessage("All done!")
}