FROM rkitchen/excerpt_base:latest
MAINTAINER Rob Kitchen <rob.kitchen@yale.edu>

##
## Add exceRpt executables
##
#ADD smallRNA_pipeline.mk /exceRpt_bin/smallRNA_pipeline
ADD exceRpt_smallRNA /exceRpt_bin/exceRpt_smallRNA
ADD mergePipelineRuns.R /exceRpt_bin/mergePipelineRuns.R
ADD mergePipelineRuns_functions.R /exceRpt_bin/mergePipelineRuns_functions.R
ADD exceRpt_Tools.jar /exceRpt_bin/exceRpt_Tools.jar
ADD LICENSE /exceRpt_bin/LICENSE
ADD README.md /exceRpt_bin/README.md

##
## Entrypoint
##
ENTRYPOINT ["make", "-f", "/exceRpt_bin/exceRpt_smallRNA", "EXE_DIR=/exceRpt_bin", "DATABASE_PATH=/exceRpt_DB", "JAVA_EXE=java", "OUTPUT_DIR=/exceRptOutput", "MAP_EXOGENOUS=off", "N_THREADS=4"]
