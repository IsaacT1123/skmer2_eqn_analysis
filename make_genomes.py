import os
import glob
import argparse

import numpy as np
from screed import fasta

def make_random_genome_with_repeats(repeat_size, contig_length, num_contigs, num_repeats, out_dir, out_name='fake_genome.fna'):
    bases = ["A", "C", "G", "T"]
    rng = np.random.default_rng()
    out_file = os.path.join(out_dir, out_name)
    print ('making random genome...')
    #remove genomes already present

    with open (out_file, 'a') as f:
        for i in range(num_contigs):
            seq = rng.choice(bases, contig_length)
            #seq = np.append(seq, np.tile(seq[-repeat_size:], num_repeats))
            f.write(f'>record{i}\n')
            f.write(f'{"".join([b for b in seq])}\n')
            print (f'made {i+1} contigs', end='\r')
    print ()
    print ('done making genome')

def mutate_fasta(in_path, out_path, d):
    reads_processed=0
    bases = set(["A", "C", "G", "T"])
    rng = np.random.default_rng()
    with open(in_path, 'r') as fin, open(out_path, 'a') as fout:
        for record in fasta.fasta_iter(fin):
            mutated_seq = "".join([rng.choice(list(bases - set([base]))) if rng.random() < d else base for base in record['sequence'].upper()])
            fout.write(f">{record['name']}_mutated\n")
            fout.write(f'{mutated_seq}\n')
            reads_processed += 1
            print (f"mutated {reads_processed} reads from {os.path.splitext(os.path.basename(in_path))[0]}...", end='\r')
    print ()
    print ("done")

def mutate_existing(args):
    for p in glob.glob(os.path.join(args.existing_dir, '*.fna')):
        out_path = os.path.join(args.out_dir, f'{os.path.splitext(os.path.basename(p))[0]}.fna')
        mutate_fasta(p, out_path, args.mutation_prob)

def mutate_random(args):
    
    for p in glob.glob(os.path.join(args.random_dir, '*')):
        os.remove(p)
    for i in range(args.num_genomes):
        make_random_genome_with_repeats(args.repeat_size, args.contig_length, args.num_contigs, args.num_repeats, args.random_dir, out_name=f'fake_genome_{i}.fna')
    for p in glob.glob(os.path.join(args.random_dir, '*.fna')):
        out_path = os.path.join(args.out_dir, f'{os.path.splitext(os.path.basename(p))[0]}.fna')
        mutate_fasta(p, out_path, args.mutation_prob)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="make some mutated genomes")

    subparsers = parser.add_subparsers()

    random_parser = subparsers.add_parser('random')
    random_parser.add_argument('--num-genomes', type=int, required=True)
    random_parser.add_argument('--repeat-size', type=int, required=True)
    random_parser.add_argument('--contig-length', type=int, required=True)
    random_parser.add_argument('--num-contigs', type=int, required=True)
    random_parser.add_argument('--num-repeats', type=int, required=True)
    random_parser.add_argument('--random-dir', type=str, required=True)
    random_parser.add_argument('--out-dir', type=str, required=True)
    random_parser.add_argument('--mutation-prob', type=float, required=True) 
    random_parser.set_defaults(func=mutate_random)

    existing_parser = subparsers.add_parser('existing')
    existing_parser.add_argument('--existing-dir', type=str, required=True)
    existing_parser.add_argument('--out-dir', type=str, required=True)
    existing_parser.add_argument('--mutation-prob', type=float, required=True) 
    existing_parser.set_defaults(func=mutate_existing)

    args = parser.parse_args()
    args.func(args)
