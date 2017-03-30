##########################################################################################
##                                                                                      ##
## Functions to combine pipeline runs for individual samples into something more useful ##
##                                                                                      ##
## Author: Rob Kitchen (r.r.kitchen@gmail.com)                                          ##
##                                                                                      ##
## Version 4.6.3 (2016-10-08)                                                           ##
##                                                                                      ##
##########################################################################################


##
## check dependencies
##
baseURL = "https://cran.r-project.org"
if(!"Rgraphviz" %in% rownames(installed.packages())) { source("http://bioconductor.org/biocLite.R"); biocLite("Rgraphviz",ask=F) }
if(!"scales" %in% rownames(installed.packages())) { install.packages("scales",repos=baseURL) }
require(Rgraphviz)
require(scales)



##
## Find the relevant sample(s) under the specified path
##
## - TODO: recurse through successive directories?
##
SearchForSampleData = function(base.dir, directory=""){
  to.return = NULL
  dir.use = paste(base.dir, directory, sep = "/")
  subdirs = dir(dir.use)

  if("ExogenousRibosomalAlignmentResults.txt" %in% subdirs)
    to.return = paste(dir.use,"ExogenousRibosomalAlignmentResults.txt",sep="")
  to.return
}



##
## Prints the given message with a timestamp
##
printMessage = function(message=""){
  cat(as.character(Sys.time()),":  ",paste(message,sep=""),"\n",sep="")
}


##
## Plots a taxonomy tree with a given set of weights
##
plotTree = function(rEG, taxonomyInfo, counts_uniq, counts_cum, title="", what){
  
  ## node parameters
  nNodes = length(nodes(rEG))
  nA <- list()
  nA$shape = rep("circle",nNodes)
  nA$fixedSize<-rep(FALSE, nNodes)
  nA$height <- nA$width <- rescale(sqrt(counts_cum/10), to=c(0.25,7))
  nA$color <- rep(rgb(0,0,0,0.25),nNodes)
  nA$style <- rep("bold", nNodes)
  if(what == "exogenousRibosomal"){
    nA$fillcolor <- sapply(counts_uniq*10, function(val){ if(val>100){val=100}; rgb(100-val,100,100-val,maxColorValue=100)})
  }else{
    nA$fillcolor <- sapply(counts_uniq*10, function(val){ if(val>100){val=100}; rgb(100-val,100-val,100,maxColorValue=100)})
  }
  
  newNodeIDs = sapply(taxonomyInfo[match(as.numeric(nodes(rEG)), taxonomyInfo$ID), ]$name, function(id){ newID=unlist(strsplit(id," ")); if(length(newID) == 1){id}else{paste(newID[1], "\n", paste(newID[-1],collapse=" "), sep="") }})
  nA$label <- paste(newNodeIDs,"\n",round(counts_cum*10)/10,"%",sep="")
  nA <- lapply(nA, function(x) { names(x) <- nodes(rEG); x})
  
  ## edge parameters
  eA <- list(arrowsize=rep(0.1,length(names(rEG@edgeData))), arrowhead=rep("none",length(names(rEG@edgeData))))
  eA <- lapply(eA, function(x) { names(x) <- names(rEG@edgeData); x})
  
  ## layout the graph
  tmp = layoutGraph(rEG, nodeAttrs=nA, edgeAttrs=eA)
  
  ## hack to make sure the node labels are visible!
  sizes = rescale(tmp@renderInfo@nodes$rWidth, to=c(0.2,1.5))
  names(sizes) = nodes(rEG)
  nodeRenderInfo(tmp) <- list(cex=sizes)
  
  graphRenderInfo(tmp) <- list(main=title)
  
  ## plot the graph
  renderGraph(tmp)
}



##
## Plot exogenous genomes
##
plotExogenousTaxonomyTrees = function(counts, cumcounts, what, output.dir, taxonomyInfo, fontScale=2, sampleGroups=NA, minPercent=0.5){
  ## add direct count to the cumulative counts matrix
  cumcounts = cumcounts+counts
  
  #counts.norm = t(t(counts*100)/colSums(counts))
  counts.norm = apply(counts, 2, function(col){ col*100/sum(col) })
  cumcounts.norm = apply(cumcounts, 2, function(col){ col*100/col[1] })
  dim(counts)
  
  ## remove nodes with < 0.1% of all reads
  #minPercent = 1
  keepRows = which(apply(counts.norm, 1, max) >= minPercent)
  keepRows = sort(unique(c(keepRows, which(apply(cumcounts.norm, 1, max) >= minPercent))))
  
  # use only paths through the tree that capture above a certain fraction of reads
  counts = counts[keepRows, , drop=F]
  cumcounts = cumcounts[keepRows, , drop=F]
  nrow(counts)
  #data_uniq = counts.norm[keepRows, , drop=F]
  #data_cum = cumcounts.norm[keepRows, , drop=F]
  #nrow(data_cum)
  
  ## Re-scale the node percentages after trimming branches to make the numbers make more sense - shouldn't make much diff to the cumcounts
  data_uniq = apply(counts, 2, function(col){ col*100/sum(col) })
  data_cum = apply(cumcounts, 2, function(col){ col*100/col[1] })
  
  ## remove edges with no useable counts (based on minPercent threshold)
  taxonomyInfo = taxonomyInfo[taxonomyInfo$ID %in% rownames(data_cum), ]
  
  ## Build the graph object
  rEG <<- new("graphNEL", nodes=as.character(taxonomyInfo$ID), edgemode="directed")
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  apply(taxonomyInfo[-1,], 1, function(row){ 
    from = trim(as.character(row[4]));
    if(from %in% taxonomyInfo$ID){ rEG <<- addEdge(trim(as.character(row[4])), trim(as.character(row[3])), rEG, 1) }
    NULL })
  
  
  data_uniq = data_uniq[match(taxonomyInfo$ID, rownames(data_uniq)), , drop=F]
  data_cum = data_cum[match(taxonomyInfo$ID, rownames(data_cum)), , drop=F]
  data_uniq[is.na(data_uniq)] = 0
  data_cum[is.na(data_cum)] = 0
  
  
  ##
  ## Write to PDF
  ##
  ## plot an average tree over all samples
  printMessage(c("Plotting a taxonomy tree based on the average of all samples "))
  pdf(file=paste(output.dir,"/exceRpt_",what,"_TaxonomyTrees_aggregateSamples.pdf",sep=""),height=7,width=15)
  plotTree(rEG, taxonomyInfo, apply(data_uniq, 1, max), rowMeans(data_cum), what=what)
  dev.off()
  
  ## plot samples individually
  printMessage(c("Plotting a separate taxonomy tree for each sample"))
  pdf(file=paste(output.dir,"/exceRpt_",what,"_TaxonomyTrees_perSample.pdf",sep=""), height=7, width=15)
  for(i in 1:ncol(data_uniq))
    plotTree(rEG, taxonomyInfo, data_uniq[,i], data_cum[,i], title=paste(colnames(data_uniq)[i]," (total reads: ",cumcounts[1,i],")", sep=""), what=what)
  dev.off()
  
  ## if there are groups of samples
  if(is.data.frame(sampleGroups)){
    printMessage(c("Plotting a separate taxonomy tree for each sample-group"))
    pdf(file=paste(output.dir,"/exceRpt_",what,"_TaxonomyTrees_perGroup.pdf",sep=""), height=7, width=15)
    for(thisgroup in levels(as.factor(sampleGroups$sampleGroup))){
      tmpDat_uniq = rowMeans(data_uniq[, match(sampleGroups[sampleGroups$sampleGroup %in% thisgroup, ]$sampleID, colnames(data_uniq)), drop=F])
      tmpDat_cum = rowMeans(data_cum[, match(sampleGroups[sampleGroups$sampleGroup %in% thisgroup, ]$sampleID, colnames(data_cum)), drop=F])
      plotTree(rEG, taxonomyInfo, tmpDat_uniq, tmpDat_cum, title=paste(thisgroup,sep=""), what=what)
    }
    dev.off()
  }
  
}


