#!/bin/bash

while getopts "d:c:o:" opts
do
	case "$opts" in
		d) kmc_db_dir="$OPTARG";;
		c) kmincov="$OPTARG";;
		o) out_dir="$OPTARG";;
		[?]) echo "invalid param"; exit 1;;
	esac	
done

if [ -z "$kmc_db_dir" ] || [ -z "$kmincov" ] || [ -z "$out_dir" ]
then
	echo "missing input param" 
	exit 1
fi



eval "$(conda shell.bash hook)"
conda activate skim-env

while read -r kmc_db
do
	db_basename=$(basename $kmc_db)
	db_name="${db_basename%.*}"
	echo "counting k-mers in db ${kmc_db_dir}/${db_name}"
	kmc_tools transform ${kmc_db_dir}/${db_name} histogram ${out_dir}/${db_name}.txt
	echo "done"
done <<<$(ls -d ${kmc_db_dir}/*.kmc_pre)

#compute kmer counts from kmer dbs

conda deactivate

