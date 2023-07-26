import os
import glob
import argparse

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import brentq, brenth

def get_genome_names(genome_dir):
    return [os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(genome_dir, "*"))]

def get_hist_orig(d, k, genome_name):
    hist_orig_path = os.path.join(args.sim_dir, "parent_child_more_complex", f"sim_more_complex_{k}_{str(d)[2:]}", "hists", "original", f"{genome_name}.hist")
    hist_orig = pd.read_csv(hist_orig_path, sep='\t', header=None).iloc[:,1] 
    return hist_orig

def get_hist_orig_skim(phred, cov, genome_name):
    hist_orig_path = os.path.join(args.genome_skim_dir, f"more_complex_phred_{phred}", f"{genome_name}_skim_{cov}x_phred_{phred}.hist")
    hist_orig = pd.read_csv(hist_orig_path, sep='\t', header=None).iloc[:,1] 
    return hist_orig

def get_hist_orig_skim_est(phred, cov, genome_name):
    hist_path = os.path.join(args.genome_respect_dir, f"{genome_name}_skim_{cov}x_phred_{phred}", "estimated-spectra.txt")
    hist_est = pd.read_csv(hist_path, sep='\t', header=0).drop(labels='sample', axis=1).transpose().iloc[:,0] 
    return hist_est

def get_params_orig_est(phred, cov, genome_name):
    hist_path = os.path.join(args.genome_respect_dir, f"{genome_name}_skim_{cov}x_phred_{phred}", "estimated-parameters.txt")
    p = pd.read_csv(hist_path, sep='\t', header=0)    
    return [p[x].squeeze() for x in ['coverage', 'sequencing_error_rate']]

def get_hist_mut(d, k, genome_name):
    hist_mut_path = os.path.join(args.sim_dir, "parent_child_more_complex", f"sim_more_complex_{k}_{str(d)[2:]}", "hists", "mutated", f"{genome_name}.txt")
    hist_mut = pd.read_csv(hist_mut_path, sep='\t', header=None).iloc[:,1] 
    return hist_mut

def get_hist_mut_skim(d, k, phred, cov, genome_name):
    hist_mut_path = os.path.join(args.sim_dir, "parent_child_more_complex", f"sim_more_complex_{k}_{str(d)[2:]}", "hists", "skims", f"mutated_phred_{phred}", f"{genome_name}_skim_{cov}x_phred_{phred}.hist")
    hist_mut = pd.read_csv(hist_mut_path, sep='\t', header=None).iloc[:,1] 
    return hist_mut

def get_hist_mut_skim_est(d, k, phred, cov, genome_name):
    hist_path = os.path.join(args.sim_dir, "parent_child_more_complex", f"sim_more_complex_{k}_{str(d)[2:]}", "respect_outputs", f"{genome_name}_skim_{cov}x_phred_{phred}", "estimated-spectra.txt")
    hist_est = pd.read_csv(hist_path, sep='\t', header=0).drop(labels='sample', axis=1).transpose().iloc[:,0] 
    return hist_est

def get_params_mut_est(d, k, phred, cov, genome_name):
    hist_mut_path = os.path.join(args.sim_dir, "parent_child_more_complex", f"sim_more_complex_{k}_{str(d)[2:]}", "respect_outputs", f"{genome_name}_skim_{cov}x_phred_{phred}", "estimated-parameters.txt")
    p = pd.read_csv(hist_mut_path, sep='\t', header=0)    
    return [p[x].squeeze() for x in ['coverage', 'sequencing_error_rate']]

def get_hist_skim_ins(d, k, phred, cov, genome_name):
    hist_skim_ins_path = os.path.join(args.sim_dir, "parent_child_more_complex", f"sim_more_complex_{k}_{str(d)[2:]}", "hists", "skims", f"ins_phred_{phred}", f"{genome_name}_skim_{cov}x_phred_{phred}_orig_cap_mut.hist")
    hist_skim_ins = pd.read_csv(hist_skim_ins_path, sep='\t', header=None).iloc[:, 1]
    return hist_skim_ins

def get_lam(cov, k, read_len):
    return cov * (read_len - k + 1) / read_len

def get_eta(cov, eps, k, read_len):
    lam = get_lam(cov, k, read_len)
    return 1 - np.exp(-lam * ((1-eps)**k))

def skmer_est(genome_name, d, phred, cov, k, read_len):
    L_orig_skim = get_hist_orig_skim(phred, cov, genome_name).sum()
    L_mut_skim = get_hist_mut_skim(d, k, phred, cov, genome_name).sum()
    truth_pc_skim = get_hist_skim_ins(d, k, phred, cov, genome_name).sum()
    eps = 10 ** (-phred/10)
    lam1 = get_lam(cov, k, read_len)
    lam2 = get_lam(cov, k, read_len)
    eta1 = get_eta(cov, eps, k, read_len)
    eta2 = get_eta(cov, eps, k, read_len)
    zeta1 = eta1 + lam1*(1- (1-eps)**k) 
    zeta2 = eta2 + lam2*(1- (1-eps)**k) 
    J = truth_pc_skim / (L_orig_skim + L_mut_skim - truth_pc_skim)
    skmer_d = 1 - ( ((zeta1 + zeta2)*J) / ((eta1 * eta2)*(1+J)) )** (1/k) 
    return skmer_d

def ins_parent_child(hist_orig, lam1, lam2, eps1, eps2, eta1, eta2, d, k, num_terms):
    left_over_orig = hist_orig[num_terms:].sum()
    hist_orig_new = hist_orig[:num_terms].copy()
    hist_orig_new[num_terms-1] += left_over_orig

    nonerr_term1 = 1 - np.power(1-eta1, 1 + np.arange(num_terms))
    nonerr_term2  = 1 - np.power((1-eta2*((1-d)**k)), 1 + np.arange(num_terms))
    nonerr_ins = np.dot(hist_orig_new, nonerr_term1*nonerr_term2)

    n1 = (1/(3*k))*lam1*k*eps1*((1-eps1)**(k-1))*(1 + np.arange(hist_orig.shape[0]))
    n21 = (1/(3*k))*(((1-d)**k)*lam2*k*eps2*((1-eps2)**(k-1)))*(1 + np.arange(hist_orig.shape[0]))
    n22 = (1/(3*k))*(k*d*(1-d)**(k-1))*(1-((1-((1-eps2)**k))**lam2))*(1 + np.arange(hist_orig.shape[0]))
    term1 = 1 - np.exp(-1*n1)
    term2= 1 - np.exp(-1*(n21 + n22))
    extra_ins = 3*k*np.dot(hist_orig, term1*term2)

    return np.dot([1, 1], [nonerr_ins, extra_ins])

def get_ins_fn(genome_name, d, phred, cov, k, read_len, num_terms):
    hist_orig = get_hist_orig(d, k, genome_name)
    lam1 = get_lam(cov, k, read_len)
    lam2= get_lam(cov, k, read_len)
    eps1 = 10 ** (-phred/10)
    eps2= 10 ** (-phred/10)
    eta1 = get_eta(cov, eps1, k, read_len)
    eta2 = get_eta(cov, eps2, k, read_len)
    def g(y):
       return ins_parent_child(hist_orig, lam1, lam2, eps1, eps2, eta1, eta2, y, k, num_terms) - np.sum(get_hist_skim_ins(d, k, phred, cov, genome_name))

    return g 

def get_ins_fn_est(genome_name, d, phred, cov, k, read_len, num_terms):
    hist_orig_skim_est = get_hist_orig_skim_est(phred, cov, genome_name)
    est_cov_orig, est_eps_orig = get_params_orig_est(phred, cov, genome_name)
    est_cov_mut, est_eps_mut = get_params_mut_est(d, k, phred, cov, genome_name)
    est_lam1 = get_lam(est_cov_orig, k, read_len)
    est_lam2= get_lam(est_cov_mut, k, read_len)
    est_eps1 = est_eps_orig
    est_eps2 = est_eps_mut
    est_eta1 = get_eta(est_cov_orig, est_eps1, k, read_len)
    est_eta2 = get_eta(est_cov_mut, est_eps2, k, read_len)
    def g(y):
       return ins_parent_child(hist_orig_skim_est, est_lam1, est_lam2, est_eps1, est_eps2, est_eta1, est_eta2, y, k, num_terms) - np.sum(get_hist_skim_ins(d, k, phred, cov, genome_name))

    return g 

def gather_data(args, use_est=False):
    genome_names = get_genome_names(args.genome_dir)
    distances = args.distances
    coverages = args.coverages
    phreds = args.phreds
    seq_errors = [10 ** (-(p/10)) for p in phreds]
    read_len = args.read_len
    k = args.ksize
    ind = pd.MultiIndex.from_product([genome_names, distances, phreds, coverages])
    df = pd.DataFrame({}, index=ind)
    #ins_parent_child(hist_orig, lam1, lam2, eps1, eps2, eta1, eta2, d, k, num_terms):
    sols = []
    for x in ind:
        try:
            if use_est:
                sol = brenth(get_ins_fn_est(x[0], x[1], x[2],  x[3], k, read_len, 5), 0, 1)
            else:
                sol = brenth(get_ins_fn(x[0], x[1], x[2],  x[3], k, read_len, 5), 0, 1)
        except:
            sol = np.nan
        sols.append(sol)
    df.loc[:, 'skmer2_est'] = sols
    #df.loc[:, 'skmer2_est'] = [brenth(get_ins_fn(x[0], x[1], x[2],  x[3], k, read_len, 5), 0, 1) for x in ind]
    dists = df.reset_index().loc[:, 'level_1']
    df.loc[:, 'skmer2_err'] =  (df.loc[:, 'skmer2_est'].to_numpy() - dists.to_numpy())/ dists.to_numpy()

    df.loc[:, 'skmer_est'] = [skmer_est(x[0], x[1], x[2],  x[3], k, read_len) for x in ind]
    df.loc[:, 'skmer_err'] =  (df.loc[:, 'skmer_est'].to_numpy() - dists.to_numpy())/ dists.to_numpy()
   
    df.index = df.index.set_names(['genome', 'distance', 'phred', 'coverage'])
    df.to_csv(os.path.join(args.out_dir, "skmer2_vs_skmer_stats_est.csv"), sep=',')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--genome-dir', type=str, required=True)
    parser.add_argument('--genome-skim-dir', type=str, required=True)
    parser.add_argument('--genome-respect-dir', type=str, required=True)
    parser.add_argument('--sim-dir', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--ksize', type=int, required=True)
    parser.add_argument('--distances', nargs='+', type=float, required=True)
    parser.add_argument('--coverages', nargs='+', type=int, required=True)
    parser.add_argument('--phreds', nargs='+', type=int, required=True)
    parser.add_argument('--read-len', type=int, required=True)
    args = parser.parse_args()
    print (f'genome_dir: {args.genome_dir}')
    print (f'distances: {args.distances}')
    print (f'coverages: {args.coverages}')
    print (f'phreds: {args.phreds}')
    print (f'read length: {args.read_len}')
    print (args)
    gather_data(args, use_est=True)
    #ins_error_by_term(args)
