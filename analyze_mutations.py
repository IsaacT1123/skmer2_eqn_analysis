import os
import glob
import argparse

import pandas as pd
from screed import fasta
import khmer


def count_mutations_to_nonunique(args, table_starting_size=1e9, num_tables=1):
    dfs = []
    for in_path in glob.glob(os.path.join(args.genome_dir_original, '*')):
        genome_name = os.path.basename(in_path)
        print (f'processing {genome_name}...')
        in_path_mutated = os.path.join(args.genome_dir_mutated, genome_name)
        ktable_orig = khmer.Counttable(args.ksize, table_starting_size, num_tables)
        records_processed = 0
        num_mutated_kmers = 0
        num_mutations_to_nonunique = 0
        #count all k-mers in original genome
        print ('counting kmers in original genome...')
        ktable_orig.consume_seqfile(in_path)
        print ('done')
        print ('counting number of mutations to nonunique kmers in mutated genome...')
        with open(in_path, 'r') as forig, open(in_path_mutated, 'r') as fmut:
            for rec_orig, rec_mut in zip(fasta.fasta_iter(forig), fasta.fasta_iter(fmut)):
                origseq = rec_orig['sequence']
                mutseq = rec_mut['sequence']
                for i in range(0, len(mutseq) - args.ksize + 1):
                    orig_kmer = origseq[i:i+args.ksize]
                    mut_kmer = mutseq[i:i+args.ksize]
                    hamming_dist = sum([a != b for a,b in zip(orig_kmer, mut_kmer)])
                    if hamming_dist > 0:
                        num_mutated_kmers += 1
                        if ktable_orig.get(mut_kmer) > 0:
                            num_mutations_to_nonunique += 1
                records_processed += 1
                #print (f'processed {records_processed} records from both genomes', end='\r')
            print ()
        df = pd.DataFrame({'genome_name':[genome_name], 'num_mutated_kmers':[num_mutated_kmers], 'num_mutations_to_nonunique':[num_mutations_to_nonunique], 'prob_mutate_to_nonunique':[num_mutations_to_nonunique/num_mutated_kmers]})
        print (df)
        dfs.append(df)
        print ('done.')
    result_df = pd.concat(dfs, axis=0)
    print (result_df)
    result_df.to_csv(os.path.join(args.out_dir, 'results.csv'), sep=',', index=False)
    print ('whole thing done')
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='none')
    parser.add_argument('--genome-dir-original', type=str, required=True)
    parser.add_argument('--genome-dir-mutated', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--ksize', type=int, required=True)

    args = parser.parse_args()
    count_mutations_to_nonunique(args)


