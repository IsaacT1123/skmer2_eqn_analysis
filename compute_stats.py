import os
import subprocess
import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binom

def aggregate_counts(hist_dir, normalize=True):
    dfs = []
    for hist_path in glob.glob(os.path.join(hist_dir, '*.txt')):
        cols = ['num_occurrences', os.path.splitext(os.path.basename(hist_path))[0]] 
        df = pd.read_csv(hist_path, sep='\t', names=cols).drop(columns=cols[0])
        #normalize counts by total number of k-mers
        df[cols[1]] = df[cols[1]] / (np.dot(1+np.arange(df.shape[0]), df[cols[1]]))**int(normalize)
        dfs.append(df)
    hists = pd.concat(dfs, axis=1)
    return hists.reindex(sorted(hists.columns), axis=1)

#def ndarange(lows, highs):
#    lows = np.array(lows)
#    highs = np.array(highs)
#    shape = (lows.shape[0], np.max(highs-lows))
#    return np.array([np.pad(np.arange(lows[i], highs[i]+1), (0, abs((highs[i]-lows[i])-shape[1]))) for i in range(lows.shape[0])])

def exp_ins_unique_pc(hists, d, k):
    return hists.sum(axis=0) * ((1-d)**k)

def exp_ins_repeats_pc(hists, d, k):
    power_terms = np.power(np.full(hists.shape[0], 1-(1-d)**k), 1+np.arange(hists.shape[0]))
    eis = hists.apply(np.dot, 0, b=1-power_terms)
    return eis

def exp_ins_repeats_pc_skim(hists, d, k, cov, cov_mut, eps, eps_mut):
    aita = 1 - np.exp(cov*(1-eps)**k)
    aita_mut = 1 - np.exp(cov_mut*(1-eps_mut)**k)
    power_terms = np.power(np.full(hists.shape[0], 1-aita*aita_mut*(1-d)**k), 1+np.arange(hists.shape[0]))
    eis = hists.apply(np.dot, 0, b=1-power_terms)

def exp_ins_repeats_ha(hists_child1, hists_child2, d, k):
    power_terms_child1 = 1 - np.power(1-(1-d/2)**k, 1+np.arange(hists_child1.shape[0]))
    power_terms_child2 = 1 - np.power(1-(1-d/2)**k, 1+np.arange(hists_child2.shape[0]))

    return ((hists_child1 + hists_child2) / 2).apply(np.dot, b=power_terms_child1*power_terms_child2)

def exp_ins_repeats_ha_skim(hists_child1, hists_child2, d, k):
    pass
    
#def expected_jaccard(ins_hists, union_hists):
#    return ins_hists.to_numpy().sum(axis=0) / union_hists.to_numpy().sum(axis=0)

def observed_sizes(hists):
    return hists.sum(axis=0)

def size_errors(exp_hist, obs_hist):
    return pd.Series((exp_hist.to_numpy() - obs_hist.to_numpy()) / obs_hist.to_numpy())

def size_errors_absolute(exp_hist, obs_hist):
    return pd.Series((exp_hist.to_numpy() - obs_hist.to_numpy()))
    

def errors_from_hists_pc(args):
    print ('computing relative error between expected and observed ins/union sizes')
    hists = aggregate_counts(args.hist_dir, normalize=False)
    ins_hists = aggregate_counts(args.ins_hist_dir, normalize=False)
    #union_hists = aggregate_counts(args.union_hist_dir, normalize=False)
    print ('computing expected/observed intersection/union sizes...')
    eis_unique = exp_ins_unique_pc(hists, args.mutation_prob, args.ksize)
    eis_repeats = exp_ins_repeats_pc(hists, args.mutation_prob, args.ksize)
    ois = observed_sizes(ins_hists)
    print ('done')
    print ('computing relative error between observed/expected union/intersection sizes...')
    error_ins_unique = size_errors(eis_unique, ois)
    error_ins_repeats = size_errors(eis_repeats, ois)
    print ('done')

    print ('computing absolute error between observed/expected union/intersection sizes...')
    error_ins_unique_abs = size_errors_absolute(eis_unique, ois)
    error_ins_repeats_abs = size_errors_absolute(eis_repeats, ois)
    print ('done.')

    print ('writing results...')
    eis_unique.to_csv(os.path.join(args.errors_dir, 'exp_ins_unique_pc.csv'), sep=",", index=True) 
    eis_repeats.to_csv(os.path.join(args.errors_dir, 'exp_ins_repeats_pc.csv'), sep=",", index=True) 

    ois.to_csv(os.path.join(args.errors_dir, 'observed_ins_sizes.csv'), sep=",", index=True) 

    error_ins_unique.to_csv(os.path.join(args.errors_dir, 'error_ins_unique.csv'), sep=",", index=True) 
    error_ins_repeats.to_csv(os.path.join(args.errors_dir, 'error_ins_repeats.csv'), sep=",", index=True) 

    error_ins_unique_abs.to_csv(os.path.join(args.errors_dir, 'error_ins_unique_abs.csv'), sep=",", index=True) 
    error_ins_repeats_abs.to_csv(os.path.join(args.errors_dir, 'error_ins_repeats_abs.csv'), sep=",", index=True) 

    print ('done')

def errors_from_hists_ha(args):
    print ('computing relative error between expected and observed ins sizes (hidden ancestor model)')
    hists_child1 = aggregate_counts(args.hist_dir_child1, normalize=False)
    hists_child2 = aggregate_counts(args.hist_dir_child2, normalize=False)
    ins_hists = aggregate_counts(args.ins_hist_dir, normalize=False)
    print ('computing expected/observed intersection sizes...')
    eis_repeats = exp_ins_repeats_ha(hists_child1, hists_child2, args.mutation_prob, args.ksize)
    ois = observed_sizes(ins_hists)
    print ('done')
    print ('computing relative error between observed/expected union/intersection sizes...')
    error_ins_repeats = size_errors(eis_repeats, ois)
    print ('done')

    print ('computing absolute error between observed/expected union/intersection sizes...')
    error_ins_repeats_abs = size_errors_absolute(eis_repeats, ois)
    print ('done.')

    print ('writing results...')
    eis_repeats.to_csv(os.path.join(args.errors_dir, 'exp_ins_repeats_pc.csv'), sep=",", index=True) 

    ois.to_csv(os.path.join(args.errors_dir, 'observed_ins_sizes.csv'), sep=",", index=True) 

    error_ins_repeats.to_csv(os.path.join(args.errors_dir, 'error_ins_repeats.csv'), sep=",", index=True) 

    error_ins_repeats_abs.to_csv(os.path.join(args.errors_dir, 'error_ins_repeats_abs.csv'), sep=",", index=True) 

    print ('done')

def rep_ratio(hists):
    return hists.apply(lambda x:1-(x.iloc[0]/np.dot(1+np.arange(x.shape[0]), x)), axis=0)

def ratios_from_hists(args):
    hists = aggregate_counts(args.hist_dir) 
    ratios = rep_ratio(hists)
    print ("1 - r1/L for all histograms:")
    print (ratios)
    print (f'mean ratio: {ratios.mean()}')
    print (f'stddev: {ratios.std()}')
    print()
    ratios.to_csv(os.path.join(args.stats_dir, 'stats.txt'))
    with open (os.path.join(args.stats_dir, 'summary.txt'), 'w') as f:
        f.write(f'mean for 1 - r1/L: {ratios.mean()}\n')
        f.write(f'stddev for 1 - r1/L: {ratios.std()}\n')

def plot_hists(args):
    cols = ['kmer_frequency', 'num_kmers']
    for p in glob.glob(os.path.join(args.hist_dir, '*')):
        plot_name = os.path.splitext(os.path.basename(p))[0]
        print (f'generating kmer freq plot for {plot_name}')
        df = pd.read_csv(p, sep='\t', names=cols)
        df['num_kmers'] = df['num_kmers'] / np.dot(1+np.arange(df.shape[0]), df['num_kmers'])
        plt.title(f'kmer freq profile for {plot_name}')
        plt.xlabel(cols[0])
        plt.ylabel(cols[1])
        plt.xlim(right=args.occ_threshold)
        plt.bar(df[cols[0]][1-args.include_unique:], df[cols[1]][1-args.include_unique:])
        plt.savefig(os.path.join(args.plot_dir, f'{plot_name}.png'))
        plt.clf()

if __name__ == '__main__':
    pass
    parser = argparse.ArgumentParser(description='kmer profile analysis functions')

    subparsers = parser.add_subparsers()
    stats_parser = subparsers.add_parser('stats')
    stats_parser.add_argument('--hist-dir', type=str, required=True)
    stats_parser.add_argument('--stats-dir', type=str, required=True)
    stats_parser.set_defaults(func=ratios_from_hists)

    errors_parser = subparsers.add_parser('errors_pc')
    errors_parser.add_argument('--hist-dir', type=str, required=True)
    errors_parser.add_argument('--ins-hist-dir', type=str, required=True)
    errors_parser.add_argument('--errors-dir', type=str, required=True)
    errors_parser.add_argument('--mutation-prob', type=float, required=True)
    errors_parser.add_argument('--ksize', type=int, required=True)
    errors_parser.set_defaults(func=errors_from_hists_pc)

    errors_parser = subparsers.add_parser('errors_ha')
    errors_parser.add_argument('--hist-dir-child1', type=str, required=True)
    errors_parser.add_argument('--hist-dir-child2', type=str, required=True)
    errors_parser.add_argument('--ins-hist-dir', type=str, required=True)
    errors_parser.add_argument('--errors-dir', type=str, required=True)
    errors_parser.add_argument('--mutation-prob', type=float, required=True)
    errors_parser.add_argument('--ksize', type=int, required=True)
    errors_parser.set_defaults(func=errors_from_hists_ha)

    plots_parser = subparsers.add_parser('plot')
    plots_parser.add_argument('--hist-dir', type=str, required=True)
    plots_parser.add_argument('--occ-threshold', type=int, required=True)
    plots_parser.add_argument('--include-unique', action='store_true', required=False)
    plots_parser.add_argument('--plot-dir', type=str, required=True)
    plots_parser.set_defaults(func=plot_hists)

    args = parser.parse_args()
    print (args)
    args.func(args)
