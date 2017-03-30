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
## Find the relative path to the script containing the required functions
##
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
#if(length(script.basename) > 0){
#  functions_master <- paste(sep="/", script.basename, "functions_master.R")
#  functions_doDEX <- paste(sep="/", script.basename, "functions_doDifferentialExpression.R")
#}else{
#  functions_master = "functions_master.R"
#  functions_doDEX = "functions_doDifferentialExpression.R"
#}
#invisible(source(functions_master))
#printMessage(paste("Sourcing",functions_doDEX,"from",script.basename)); invisible(source(functions_doDEX)); cat("\n")




##
## Check inputs and do the plots
##
args<-commandArgs(TRUE)
if(length(args) == 0){
  ## if no data directory is specified, throw an error message
  cat("\nERROR: no input data directory specified!\n\n")
  cat("Usage: Rscript plotExogenous.R <data path> [output path]\n\n")
}else{
  data.dir = args[1]
  if(length(args) >= 2){
    output.dir = args[2]
  }else{
    output.dir = data.dir
  }

  paths = unique(SearchForSampleData(data.dir))
  
  if(!is.null(paths)){
    library(readr)
    data = as.data.frame(read_tsv(paths[1]))
    
    counts = data.frame(s1=data$readCount_direct)
    cumcounts = data.frame(s1=data$readCount_inherited)
    rownames(counts) = rownames(cumcounts) = data$ID
    
    plotExogenousTaxonomyTrees(counts, cumcounts, what="s1", output.dir, taxonomyInfo=data[,1:5], fontScale=2, sampleGroups=NA, minPercent=0.5)
  }
}



