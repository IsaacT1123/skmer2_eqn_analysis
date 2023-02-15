#!/bin/bash

while getopts "g:s:k:p:t:" opts 
do
	case "$opts" in
		g) genome_dir="$OPTARG";;
		s) out_dir="$OPTARG";;
		k) ksize="$OPTARG";;
		p) mutation_prob="$OPTARG";;
		t) threads="$OPTARG";;
		[?]) echo "invalid input param"; exit 1;;
	esac
done


if [ -z "$genome_dir" ] || [ -z "$out_dir" ] || [ -z "$ksize" ] || [ -z "$mutation_prob" ] || [ -z "$threads" ]
then
	echo "missing input parameter"
	exit 1
fi

scripts_dir=$(dirname $(realpath "$0") )
echo "scripts directory is $scripts_dir"

genome_dir_basename=$(basename $genome_dir)
prob_cut=${mutation_prob##*.}
sim_dir=${out_dir}/sim_fruit_flies_${ksize}_${prob_cut}
mkdir -p ${sim_dir}/results
cd $sim_dir

#mutate source genome
#get kmers from source genome
#for each kmer in source genome, compute all kmers with hamming distance of 1
#count occurrences of those kmers in mutated genome
#echo "Mutate genomes at $genome_dir with mutation probability of $mutation_prob. Using kmer length of $ksize, keeping kmers occuring at least $kmincov times, using $threads threads. Output results to $sim_dir"


eval "$(conda shell.bash hook)"
#conda activate skim-env
#
#python ${scripts_dir}/make_genomes.py existing --existing-dir $genome_dir --out-dir $sim_dir/mutated_genomes --mutation-prob $mutation_prob 
#
#conda deactivate
#
conda activate khmer-env

#estimate probability of a mutation resulting in a non-unique kmer
#python ${scripts_dir}/analyze_mutations.py --genome-dir-original $genome_dir --genome-dir-mutated ${sim_dir}/mutated_genomes --out-dir ${sim_dir}/results --ksize $ksize
python ${scripts_dir}/count_adj.py --genome-dir-original $genome_dir --out-dir ${sim_dir}/results --ksize $ksize

conda deactivate

echo "whole pipe done"
