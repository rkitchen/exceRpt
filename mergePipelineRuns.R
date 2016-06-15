#######################################################################################
##                                                                                   ##
## Script to combine pipeline runs for individual samples into something more useful ##
##                                                                                   ##
## Author: Rob Kitchen (rob.kitchen@yale.edu)                                        ##
##                                                                                   ##
## Version 3.2.0 (2015-11-02)                                                        ##
##                                                                                   ##
#######################################################################################


##
## Check inputs
##
args<-commandArgs(TRUE)
if(length(args) == 0){
  
  ## if no data directory is specified, throw an error message
  cat("\nERROR: no input data directory specified!\n\n")
  cat("Usage: Rscript mergePipelineRuns.R <data path> [output path]\n\n")

}else{
  
  data.dir = args[1]
  if(length(args) >= 2){
    output.dir = args[2]
    if(length(args) == 3){
      classifier.path = args[3]
    }
  }else{
    output.dir = data.dir
  }
  
  
  ##
  ## Find the relative path to the script containing the required functions
  ##
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  if(length(script.basename) > 0){
	  other.name <- paste(sep="/", script.basename, "mergePipelineRuns_functions.R")
  }else{
	  other.name = "mergePipelineRuns_functions.R"
  }
  print(paste("Sourcing",other.name,"from",script.name))
  source(other.name)
  cat("\n")
  
  
  ##
  ## Process all samples under this directory
  ##
  processSamplesInDir(data.dir, output.dir, scriptDir=script.basename)
}


