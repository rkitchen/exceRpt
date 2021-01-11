conda activate snakemake

##
## Run on AWS
##
PATH_ANN=/scratch/Annotations/exceRpt
#aws s3 sync s3://kitchen-mgh-public/exceRpt/DATABASE/v5.0/exceRptDB_mm10 $PATH_ANN/
aws s3 sync s3://kitchen-mgh-public/exceRpt/DATABASE/v5.0/exceRptDB_hg38 $PATH_ANN/

BASE=/efs/Pipelines/exceRpt
#PATH_OUT=$BASE/TESTING/output
PATH_OUT=/scratch/exceRpt_runs
mkdir -p $PATH_OUT

conda activate snakemake

## Visualise pipeline
snakemake -s $BASE/exceRpt_smallRNA.py --cores 8 \
	--directory=$PATH_OUT --configfile="$BASE/TESTING/config_aws.yaml" --dag \
	| dot -Tpdf > $BASE/TESTING/dag_full.pdf
snakemake -s $BASE/exceRpt_smallRNA.py --cores 8 \
	--directory=$PATH_OUT --configfile="$BASE/TESTING/config_aws.yaml" --rulegraph \
	| dot -Tpdf > $BASE/TESTING/dag_simple.pdf

## Run pipeline
snakemake -s $BASE/exceRpt_smallRNA.py --cores 2 \
	--directory=$PATH_OUT --configfile="$BASE/TESTING/config_aws.yaml" -n -p

## Tidy up
rm -r $PATH_OUT
