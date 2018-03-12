exceRpt small-RNAseq pipeline
====================

Software for preprocessing, filtering, alignment, and reporting of smallRNA-seq datasets.


AUTHOR/SUPPORT:

Rob Kitchen, r [dot] r [dot] kitchen [@] gmail



CONTENTS:

exceRpt_smallRNA -- This choreographs the processing, filtering, and alignment of a single smallRNA-seq sample. The script is a makefile and 

mergePipelineRuns.R -- This script will take as input a directory containing 1 or more subdirectories or zipfiles containing output from the pipeline above. In this way, results from 1 or more smallRNA-seq samples can be combined, several QC plots are generated, and the read-counts are normalised ready for downstream analysis by clustering and/or differential expression.

Please see the exceRpt mainpage (https://rkitchen.github.io/exceRpt) for instructions as to how to use the software


INSTALLATION:

exceRpt_smallRNA -- requires numerous dependences that require some familiarity with UNIX.  The software is installed on the Genboree Workbench (www.genboree.org) which provides a graphical interface to the pipeline that is free for academic use.  Alternatively, development of a Docker image (https://hub.docker.com/r/rkitchen/excerpt) for this software is available, if you are interested in this please feel free browse the documentation (https://rkitchen.github.io/exceRpt) or get in touch using the email address above.

mergePipelineRuns.R -- is comparatively much simpler to install.  Once the R software (http://cran.r-project.org/) is set up on your system the script should automatically identify and install all required dependencies.  Again, this script is available on the Genboree Workbench (www.genboree.org) and is also free for academic use.
