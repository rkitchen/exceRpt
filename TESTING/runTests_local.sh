conda activate snakemake

##
## Run locally
##
BASE=~/WORK/exceRpt/exceRpt_github
PATH_OUT=$BASE/TESTING/output
mkdir -p $PATH_OUT

## Visualise pipeline
snakemake -s $BASE/exceRpt_smallRNA.py --cores 8 \
	--directory=$PATH_OUT --configfile="$BASE/TESTING/config_local.yaml" --dag \
	| dot -Tpdf > $BASE/TESTING/dag_full.pdf
snakemake -s $BASE/exceRpt_smallRNA.py --cores 8 \
	--directory=$PATH_OUT --configfile="$BASE/TESTING/config_local.yaml" --rulegraph \
	| dot -Tpdf > $BASE/TESTING/dag_simple.pdf

## Run pipeline
snakemake -s $BASE/exceRpt_smallRNA.py --cores 2 \
	--directory=$PATH_OUT --configfile="$BASE/TESTING/config_local.yaml" -n -p

## Tidy up
rm -r $PATH_OUT