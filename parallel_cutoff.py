import numpy as np
from math import *
import time
from utils import *

from joblib import Parallel, delayed
# what are your inputs, and what operation do you want to
# perform on each input. For example...
num_cores = 12#multiprocessing.cpu_count()                                     

n = 4 # number of segments (i.e., no. of repeaters -1 )
Nm = 1e6
F_link = 1 # fidelity of local BP
mu_link = 1 # depolarizing noise channel parameter (1: no noise, 0: fully depolarized)

τ_coh_list = np.logspace(-4.3,-2,10) # coherence time [sec]
Le2e_list = np.linspace(50,400,4)

num_τ = 10

def runner(i_L):

    raw_rate_par = np.zeros((len(τ_coh_list),num_τ))
    skr_par = np.zeros((len(τ_coh_list),num_τ))
    Fe2e_par = np.zeros((len(τ_coh_list),num_τ))

    skr_par_opt = np.zeros((len(τ_coh_list)))
    raw_rate_par_opt = np.zeros((len(τ_coh_list)))
    Fe2e_par_opt = np.zeros(len(τ_coh_list))
    τ_cut_opt =  np.zeros(len(τ_coh_list))
    skr_par_no_cut =  np.zeros(len(τ_coh_list))
    raw_rate_par_no_cut = np.zeros(len(τ_coh_list))
    Fe2e_par_no_cut = np.zeros(len(τ_coh_list))

    tic = time.time()
    Le2e = Le2e_list[i_L]
    Ls = [Le2e/n]*n
    τ_cut_list = np.logspace(-0.5,2,num_τ)*Ls[0]/c/2 # cutoff [sec]

    for i_coh, τ_coh in enumerate(τ_coh_list):
        for i_t, τ_cut in enumerate(τ_cut_list):
            raw_rate_par[i_coh,i_t], skr_par[i_coh,i_t], Fe2e_par[i_coh,i_t] = T_parallel_cutoff(τ_cut, τ_coh, mu_link, F_link, Ls, cct= True, Nmax=Nm)
            if isnan(skr_par[i_coh,i_t]):
                skr_par[i_coh,i_t] = 0#1e-20

        idx = np.argmax(skr_par[i_coh,:])
        skr_par_opt[i_coh] = skr_par[i_coh,idx]
        raw_rate_par_opt[i_coh] = raw_rate_par[i_coh,idx]
        Fe2e_par_opt[i_coh] = Fe2e_par[i_coh,idx]
        
        τ_cut_opt[i_coh] = τ_cut_list[idx]
        raw_rate_par_no_cut[i_coh], skr_par_no_cut[i_coh], Fe2e_par_no_cut[i_coh] = T_parallel_no_cutoff(τ_coh, mu_link, F_link, Ls, cct= True, Nmax=Nm)

    out_dir = 'data_color_plot/par/'
    fname = f"n_{n}_L_{i_L}.npz"
    np.savez(out_dir+fname, τ_coh_list, Le2e_list,
            raw_rate_par, skr_par, Fe2e_par, 
            skr_par_opt, raw_rate_par_opt, Fe2e_par_opt, τ_cut_opt,
            raw_rate_par_no_cut, skr_par_no_cut, Fe2e_par_no_cut)
    toc = time.time()
    print(i_L, ", Elapsed time: ", toc-tic)
    return 0


results = Parallel(n_jobs=num_cores)(delayed(runner)(i_L) for i_L in range(len(Le2e_list)))
