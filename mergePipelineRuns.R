#######################################################################################
##                                                                                   ##
## Script to combine pipeline runs for individual samples into something more useful ##
##                                                                                   ##
## Author: Rob Kitchen (rob.kitchen@yale.edu)                                        ##
##                                                                                   ##
## Version 3.1.1 (2015-11-01)                                                        ##
##                                                                                   ##
#######################################################################################


##
## Define inputs:
##
args<-commandArgs(TRUE)

##
## Check inputs
##
#stopifnot(length(args) >= 1)
if(length(args) >= 1){
  data.dir = args[1]
  if(length(args) >= 2){
    output.dir = args[2]
    if(length(args) == 3){
      classifier.path = args[3]
    }else{
      #classifier.path = "/Users/robk/Box Sync/Work/miRNA/DataSets"
      classifier.path = "/Users/robk/Box Sync/Work/exRNA/Pipeline/smallRNA_postprocessing_code/tissueClassifier.RData"      
    }
  }else{
    output.dir = data.dir
  }
}else{ 
  ## if no data directory is specified, throw an error message and try a hard coded path (for testing)
  cat("ERROR: no input data directory specified!\n\n")
  cat("Usage: Rscript mergePipelineRuns.R <data path> [output path]\n\n")
  #data.dir = "~/Box Sync/Work for other people/FraminghamSmallRNA/Data"
  #output.dir = "~/Box Sync/Work for other people/FraminghamSmallRNA"
  #data.dir = "~/WORK/YALE_offline/miRNA/PipelineOutput"
  #data.dir = "~/WORK/YALE_offline/exRNA/newPipelineTest"  
  #data.dir = "~/WORK/YALE_offline/exRNA/DavidGalas"
  #data.dir = "~/WORK/YALE_offline/exRNA/TomTuschl"
  #data.dir = "~/WORK/YALE_offline/exRNA/TomTuschl_singleSample"
  #data.dir = "~/WORK/YALE_offline/exRNA/AlCharest"	
  #data.dir = "~/WORK/YALE_offline/exRNA/LouiseLaurent"
  #data.dir = "~/WORK/YALE_offline/exRNA/YanAssman"
  #data.dir = "~/WORK/YALE_offline/exRNA/exceRptLibraryAnalysis"
  #data.dir = "~/WORK/YALE_offline/exRNA/exceRptCombinatorialAnalysis"
  #data.dir = "~/WORK/YALE_offline/exRNA/TESTING"
  #data.dir = "~/WORK/YALE_offline/exRNA/TESTING/standardTests"
  #data.dir = "~/Downloads/JAMES_TEST"
  #data.dir = "~/WORK/YALE_offline/exRNA/AmyBuck/Results"
  #data.dir = "~/WORK/YALE_offline/BrainSpan/smallRNA_new/exceRpt_Results"
  output.dir = data.dir
}



source('~/Box Sync/Work/exRNA/Pipeline/smallRNApipe/mergePipelineRuns_functions.R')
processSamplesInDir(data.dir, output.dir)
