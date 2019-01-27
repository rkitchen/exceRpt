FROM rkitchen/excerpt_base:latest
MAINTAINER Rob Kitchen <r.r.kitchen@gmail.com>

##
## Add exceRpt executables
##
ADD exceRpt_smallRNA /exceRpt_bin/exceRpt_smallRNA
ADD mergePipelineRuns.R /exceRpt_bin/mergePipelineRuns.R
ADD mergePipelineRuns_functions.R /exceRpt_bin/mergePipelineRuns_functions.R
ADD exceRpt_Tools.jar /exceRpt_bin/exceRpt_Tools.jar
ADD LICENSE /exceRpt_bin/LICENSE
ADD README.md /exceRpt_bin/README.md

##
## Add baseDB and example raw data
##
ADD exceRpt_coreDB /exceRpt_DB/
ADD ExampleData/testData_human.fastq.gz /exceRptInput/testData_human.fastq.gz


##
## Entrypoint
##
ENTRYPOINT ["make", "-f", "/exceRpt_bin/exceRpt_smallRNA", "EXE_DIR=/exceRpt_bin", "DATABASE_PATH=/exceRpt_DB", "JAVA_EXE=java", "OUTPUT_DIR=/exceRptOutput", "MAP_EXOGENOUS=off", "N_THREADS=4"]


