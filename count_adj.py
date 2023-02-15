import os
import glob
import argparse

import pandas as pd
from screed import fasta
import khmer

def get_one_adj(kmer):
    kmers = []
    alphabet = set(['A', 'C', 'G', 'T'])
    bases = list(kmer)
    for i in range(len(kmer)):
        for c in alphabet - set(kmer[i]):
            bases[i] = c
            kmers.append("".join(bases))
            bases[i] = kmer[i]
    return set(kmers)
    

def prob_mutation_is_nonunique(args, table_starting_size=1e9, num_tables=1):
    dfs = []
    for in_path in glob.glob(os.path.join(args.genome_dir_original, '*')):
        genome_name = os.path.basename(in_path)
        seen_orig = set()
        print (f'processing {genome_name}...')
        print ('getting kmers from original genome...')
        with open(in_path, 'r') as forig:
            for rec_orig in fasta.fasta_iter(forig):
                origseq = rec_orig['sequence']
                for i in range(0, len(origseq) - args.ksize + 1):
                    curr_kmer = origseq[i:i+args.ksize] 
                    seen_orig.add(curr_kmer)
        print ('done')
        running_prob = 0
        print ('counting adjacent kmers also in genome...')
        for kmer in seen_orig:
            adj = get_one_adj(kmer)
            adj_in_orig = adj.intersection(seen_orig)
            running_prob += len(adj_in_orig) / len(adj)
        print ('done.') 
        df = pd.DataFrame({'genome_name':[genome_name], 'avg_prob_mutate_to_nonunique':[running_prob / len(seen_orig)]})
        print (df)
        dfs.append(df)
        print ('done this iteration.')
        break
    
    result_df = pd.concat(dfs, axis=0)
    print (result_df)
    result_df.to_csv(os.path.join(args.out_dir, 'prob_estimate_results.csv'), sep=',', index=False)
    print ('whole thing done')
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='none')
    parser.add_argument('--genome-dir-original', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--ksize', type=int, required=True)

    args = parser.parse_args()
    prob_mutation_is_nonunique(args)


