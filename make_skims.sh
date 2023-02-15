#!/bin/bash
set -x 
while getopts "g:s:c:p:l:" opts 
do
	case "$opts" in
		g) genome_dir="$OPTARG";;
		s) skim_dir="$OPTARG";;
		c) coverage="$OPTARG";;
		p) phred_quality="$OPTARG";;
		l) read_len="$OPTARG";;
		[?]) echo "invalid input param"; exit 1;;
	esac
done
for p in $(ls -d ${genome_dir}/*)
do
    pbasename=$(basename $p)
    pfilename=${pbasename%%.*}
    pext=${pbasename##*.}
    echo "generating skim of $pfilename at ${coverage}x coverage with phred quality of ${phred_quality}..."
    /calab_data/mirarab/home/isthomas/dev/art_bin_MountRainier/art_illumina -q -i $p -l $read_len -f $coverage -qL ${phred_quality} -qU ${phred_quality} -o ${pfilename}_skim_${coverage}x
    echo "done"
done
