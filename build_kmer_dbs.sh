#!/bin/bash

while getopts "g:k:c:j:o:t:" opts
do
	case "$opts" in
		g) genomes_dir="$OPTARG";;
		k) ksize="$OPTARG";;
		c) kmincov="$OPTARG";;
		j) tmp_dir="$OPTARG";;
		o) out_dir="$OPTARG";;
		t) threads="$OPTARG";;
		[?]) echo "invalid param"; exit;;
	esac
done


if [ -z "$genomes_dir" ] || [ -z "$ksize" ] || [ -z "$kmincov" ] || [ -z "$tmp_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ]	
then
	echo "missing input param" 
	exit 1
fi

eval "$(conda shell.bash hook)"
conda activate skim-env

while read -r genome_path
do
	genome_basename=$(basename $genome_path)
	genome_name="${genome_basename%.*}"
	echo "generating kmer db for $genome_name"
	kmc -k$ksize -ci$kmincov -t$threads -fm $genome_path ${out_dir}/${genome_name} ${tmp_dir}
	echo "done"
done <<<$(ls -d ${genomes_dir}/*)  

conda deactivate
