#!/bin/bash

while getopts "g:b:s:k:c:p:x:t:" opts 
do
	case "$opts" in
		g) genome_dir="$OPTARG";;
		b) original_db_dir="$OPTARG";;
		s) out_dir="$OPTARG";;
		k) ksize="$OPTARG";;
		c) kmincov="$OPTARG";;
		p) mutation_prob="$OPTARG";;
		x) occ_threshold="$OPTARG";;
		t) threads="$OPTARG";;
		[?]) echo "invalid input param"; exit 1;;
	esac
done


if [ -z "$genome_dir" ] || [ -z "$original_db_dir" ] || [ -z "$out_dir" ] || [ -z "$ksize" ] || [ -z "$kmincov" ] || [ -z "$mutation_prob" ] || [ -z "$occ_threshold" ] || [ -z "$threads" ]
then
	echo "missing input parameter"
	exit 1
fi

scripts_dir=$(dirname $(realpath "$0") )
echo "scripts directory is $scripts_dir"
echo "original db dir is $original_db_dir"

genome_dir_basename=$(basename $genome_dir)
prob_cut=${mutation_prob##*.}
sim_dir=${out_dir}/sim_${genome_dir_basename}_${ksize}_${prob_cut}
mkdir -p ${sim_dir}/tmp/{kmc,respect,mash} ${sim_dir}/mutated_genomes/{child1,child2} ${sim_dir}/{dbs,hists,stats,plots}/{original,ins,union} ${sim_dir}/{dbs,hists,stats,plots}/mutated/{child1,child2} ${sim_dir}/respect_output/original ${sim_dir}/respect_output/mutated/{child1,child2}
mkdir -p ${sim_dir}/skmer_pop_input/{genomes,params} ${sim_dir}/skmer_output
mkdir -p ${sim_dir}/skmer_pop_output
mkdir -p ${sim_dir}/mash_output
mkdir -p ${sim_dir}/errors
ls ${sim_dir}
cd $sim_dir

echo "process  genomes in $genome_dir with mutation probability of $mutation_prob. Using kmer length of $ksize, keeping kmers occuring at least $kmincov times, using $threads threads. Output results to $sim_dir"


eval "$(conda shell.bash hook)"
conda activate skim-env

echo "mutating input genomes into child genomes half the original distance away with original mutation probability of $mutation_prob"
#python ${scripts_dir}/make_genomes.py random --num-genomes 1 --repeat-size 0 --contig-length 100000  --num-contigs 10000 --num-repeats 0 --random-dir $genome_dir --out-dir $sim_dir/mutated_genomes --mutation-prob $mutation_prob 
divisor=2
dist_half=$(echo "scale=4 ; $mutation_prob / $divisor" | bc)
echo "making child1..."
python ${scripts_dir}/make_genomes.py existing --existing-dir $genome_dir --out-dir ${sim_dir}/mutated_genomes/child1 --mutation-prob $dist_half 
echo "done making child1."

echo "making child2..."
python ${scripts_dir}/make_genomes.py existing --existing-dir $genome_dir --out-dir ${sim_dir}/mutated_genomes/child2 --mutation-prob $dist_half 
echo "done making child2."

#echo "computing mash jaccard"
#for p in $(ls -d ${genome_dir}/*)
#do
#	p_basename=$(basename $p)
#	p_name="${p_basename%.*}"
#	mash sketch -k $ksize -s 100000 -o ${sim_dir}/tmp/mash/${p_name} $p
#	mash sketch -k $ksize -s 100000  -o ${sim_dir}/tmp/mash/${p_name}_mutated ${sim_dir}/mutated_genomes/${p_name}_mutated.fna 
#	echo "computing distance..."
#	mash dist ${sim_dir}/tmp/mash/${p_name}.msh ${sim_dir}/tmp/mash/${p_name}_mutated.msh > ${sim_dir}/mash_output/dist_${p_name}.txt 
#	echo "done"
#
#done

conda deactivate

echo "building kmer dbs from mutated genomes"
bash ${scripts_dir}/build_kmer_dbs.sh -g ${sim_dir}/mutated_genomes/child1 -k $ksize -c $kmincov -j ${sim_dir}/tmp/kmc -o ${sim_dir}/dbs/mutated/child1 -t $threads
bash ${scripts_dir}/build_kmer_dbs.sh -g ${sim_dir}/mutated_genomes/child2 -k $ksize -c $kmincov -j ${sim_dir}/tmp/kmc -o ${sim_dir}/dbs/mutated/child2 -t $threads
echo "done"

echo "intersecting dbs"
bash ${scripts_dir}/combine_dbs_ha.sh -a ${sim_dir}/dbs/mutated/child1 -b ${sim_dir}/dbs/mutated/child2 -c $kmincov -e intersect -o ${sim_dir}/dbs/ins 
echo "done"

echo "generating histogram data for all kmc dbs"
bash ${scripts_dir}/count_kmers.sh -d $original_db_dir -c $kmincov -o ${sim_dir}/hists/original 
bash ${scripts_dir}/count_kmers.sh -d ${sim_dir}/dbs/mutated/child1 -c $kmincov -o ${sim_dir}/hists/mutated/child1
bash ${scripts_dir}/count_kmers.sh -d ${sim_dir}/dbs/mutated/child2 -c $kmincov -o ${sim_dir}/hists/mutated/child2
bash ${scripts_dir}/count_kmers.sh -d ${sim_dir}/dbs/ins -c $kmincov -o ${sim_dir}/hists/ins
echo "done"

conda activate skim-env

echo "computing stats for kmer histograms"
python ${scripts_dir}/compute_stats.py stats --hist-dir ${sim_dir}/hists/original --stats-dir ${sim_dir}/stats/original
python ${scripts_dir}/compute_stats.py stats --hist-dir ${sim_dir}/hists/mutated/child1 --stats-dir ${sim_dir}/stats/mutated/child1
python ${scripts_dir}/compute_stats.py stats --hist-dir ${sim_dir}/hists/mutated/child2 --stats-dir ${sim_dir}/stats/mutated/child2
python ${scripts_dir}/compute_stats.py stats --hist-dir ${sim_dir}/hists/ins --stats-dir ${sim_dir}/stats/ins
echo "done"

echo "computing error between expected & actual intersection sizes"
python ${scripts_dir}/compute_stats.py errors_ha --hist-dir-child1 ${sim_dir}/hists/mutated/child1 --hist-dir-child2 ${sim_dir}/hists/mutated/child2 --ins-hist-dir ${sim_dir}/hists/ins --mutation-prob $mutation_prob --ksize $ksize --errors-dir ${sim_dir}/errors
echo "done"

cd $scripts_dir

#rm ${sim_dir}/skmer_output/mutated/*
#rm ${sim_dir}/skmer_output/original/*
#rm ${sim_dir}/skmer_pop_output/mutated/*
#rm ${sim_dir}/skmer_pop_output/original/*
#for p in $(ls -d ${genome_dir}/*)
#do
#	rm ${sim_dir}/skmer_pop_input/genomes/*
#	rm ${sim_dir}/skmer_pop_input/params/*
#	p_basename=$(basename $p)
#	p_name="${p_basename%.*}"
#	#copy genome and mutated counterpart to skmer input genomes directory
#	#estimate spectra and params with respect and copy to skmer input params directory
#	echo "computing respect spectra & params and estimating skmer/skmer-pop distance between ${p_name} and ${p_name}_mutated..."
#	cp ${genome_dir}/${p_name}.fna ${sim_dir}/skmer_pop_input/genomes/original.fna 
#	cp ${sim_dir}/mutated_genomes/${p_name}_mutated.fna ${sim_dir}/skmer_pop_input/genomes/mutated.fna
#	set -x
#	echo "running skmer..."
#	skmer reference -l ${sim_dir}/skmer_output -o ${sim_dir}/skmer_output/ref-dist-mat-${p_name} ${sim_dir}/skmer_pop_input/genomes 
#	echo "done."
#	echo "running respect and skmer-pop..."
#	respect -d ${sim_dir}/skmer_pop_input/genomes -N 1000 --tmp ${sim_dir}/tmp/respect -o ${sim_dir}/skmer_pop_input/params
#	set +x
#
#	python ${scripts_dir}/skmer-pop.py reference -l ${sim_dir}/skmer_pop_output --param ${sim_dir}/skmer_pop_input/params/estimated-parameters.txt --spec ${sim_dir}/skmer_pop_input/params/estimated-spectra.txt -o ${sim_dir}/skmer_pop_output/ref-dist-mat-${p_name} -k $ksize -p $threads ${sim_dir}/skmer_pop_input/genomes
#	echo "done"
#done

#echo "whole pipe done for this round."

conda deactivate 
