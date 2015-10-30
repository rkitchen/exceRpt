exceRpt small-RNAseq pipeline
====================

Software for preprocessing, filtering, alignment, and reporting of smallRNA-seq datasets.


AUTHOR/SUPPORT:

Rob Kitchen, rob.kitchen@yale.edu



CONTENTS:

smallRNA_pipeline -- This choreographs the processing, filtering, and alignment of a single smallRNA-seq sample. The script is a makefile and 

mergePipelineRuns -- This script will take as input a directory containing 1 or more subdirectories or zipfiles containing output from the pipeline above. In this way, results from 1 or more smallRNA-seq samples can be combined, several QC plots are generated, and the read-counts are normalised ready for downstream analysis by clustering and/or differential expression.



INSTALLATION:

smallRNA_pipeline -- requires numerous dependences that require some familiarity with UNIX.  The software is installed on the Genboree Workbench (www.genboree.org) which provides a graphical interface to the pipeline that is free for academic use.  Alternatively, development of a Docker image for this software is ongoing, if you are interested in this please feel free to get in touch using the email address above.

mergePipelineRuns -- is comparatively much simpler to install.  Once the R software (http://cran.r-project.org/) is set up on your system the script should automatically identify and install all required dependencies.  Again, this script is available on the Genboree Workbench (www.genboree.org) and is also free for academic use.
