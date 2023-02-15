#!/bin/bash

while getopts "a:b:c:e:o:" opts
do
	case "$opts" in
		a) db_dir1="$OPTARG";;
		b) db_dir2="$OPTARG";;
		c) kmincov="$OPTARG";;
		e) operation="$OPTARG";;
		o) out_dir="$OPTARG";;
		[?]) echo "invalid param"; exit;;
	esac
done


if [ -z "$db_dir1" ] || [ -z "$db_dir2" ] || [ -z "$kmincov" ] || [ -z "$operation" ] || [ -z "$out_dir" ]
then
	echo "missing input param" 
	exit 1
fi

eval "$(conda shell.bash hook)"
conda activate skim-env

for p1 in $(ls -d ${db_dir1}/*.kmc_pre)
do
	p1_basename=$(basename $p1)
	p1_name="${p1_basename%.*}"
	echo "peforming $operation on dbs from $p1_name and ${p1_name}_mutated"
	kmc_tools simple ${db_dir1}/${p1_name} -ci$kmincov ${db_dir2}/${p1_name} -ci$kmincov $operation ${out_dir}/${p1_name}_child1_cap_child2
	echo "done"
done

conda deactivate
