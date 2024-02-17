import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from math import *

# geometric distribution for trials prob(n) = p q^(n-1) where p is success probability
F_geo = lambda x,p: np.floor(np.log(1-x)/np.log(1-p))
f = lambda x,p: p* (1-p)**(x-1)
# photon transmission probability in fiber (i.e., 0.2dB/km)
Trans = lambda x: 10**(-0.2*x/10)
# Binary Shanon entropy
def h(p_list):
    y_list = np.zeros(len(p_list))
    for i,p in enumerate(p_list):
        if p<1e-6 or (1-p)<1e-6:
            y_list[i]= 0
        else:
            y_list[i]= -p*np.log2(p)-(1-p)*np.log2(1-p)
    return y_list

# constants
c = 2e5 # speed of light in fiber [km/s]
p_link = 1.0 # photon insertion loss incorporates various efficiencies of the experimental hardware

def T_parallel_no_cutoff(τ_coh, mu_link, F_link, links, cct= True, Nmax=10000):
    """ Runs Monte-Carlo iterations to calculate performance metrics for asynchronous parallel scheme
    inputs:
        τ_coh: coherence time of quantum memories
        mu_link: parameter in 2qubit depolarizing channel describing noisy link-level entanglement and
        entanglement swapping error
        F_link: fidelity of link level entanglement (i.e.,quality of locally generated Bell pairs)
        links: list of segment (link) lengths in km
        cct: bool flag to turn on/off classical comm.
        Nmax: No. of MC iterations (if you see oscillations in secret key curves, you need to try larger numbers)

    outputs:
        Raw_rate: 1/ expected value of total time for e2e entanglement delivery
        *** application specific quantities:
        skr: secret key rate for qkd (does not include idle times of end memories)
        F_e2e: e2e entanglement fidelity for entanglement distrubtion (does include idle times of end memories)
    """
    Nmax = int(Nmax) # to make sure Nmax is an integer
    N_links = len(links) # number of links, i.e. no. of repeaters + 1
    if type(links) != np.ndarray:
        links = np.array(links)
    τs = links/c
    Ns = np.zeros((N_links,Nmax))
    # Ns[0,:] = (2*F_geo(np.random.rand(Nmax),p_link*Trans(links[0]))-1)*τs[0]
    p_i = p_link*Trans(links[0])
    Ns[0,:] = (2*F_geo(np.random.uniform(low=p_i, high=1, size=(Nmax,)),p_i)-1)*τs[0]
    for i in range(1,N_links):
        # Ns[i,:] = 2*F_geo(np.random.rand(Nmax),p_link*Trans(links[i]))*τs[i]
        p_i = p_link*Trans(links[i])
        Ns[i,:] = 2*F_geo(np.random.uniform(low=p_i, high=1, size=(Nmax,)),p_i)*τs[i]

    Ts = np.zeros((N_links-1,Nmax))
    Ts[0,:] = np.abs(Ns[1,:]-Ns[0,:])+ 2*τs[1]
    for i in range(1,N_links-1):
        Ts[i,:] = np.abs(Ns[i+1,:]-Ns[i,:]+τs[i])+ 2*τs[i+1]
    f_memory_qkd = np.mean(np.exp(-np.sum(Ts,axis=0)/τ_coh))
    if cct: # inclduing classical comm. 
        Tsw = np.zeros((N_links-1,Nmax))
        Tsw[0,:] = np.max([Ns[1,:],Ns[0,:]],axis=0) + np.sum(τs[:1])
        for i in range(1,N_links-1):
            Tsw[i,:] = np.max([Ns[i+1,:],Ns[i,:]-τs[i]],axis=0)+ np.sum(τs[:(i+1)])
            
        T_tot = np.mean(np.max(Tsw,axis=0))           
        Ts_A = np.max(Tsw,axis=0)-Ns[0,:] + τs[0] 
        Ts_B = np.max(Tsw,axis=0)-Ns[-1,:] + τs[-1] 
        f_memory_bell = np.mean(np.exp(-(np.sum(Ts,axis=0)+Ts_A+Ts_B)/τ_coh))
    else:
        T_tot = np.mean(np.max(Ns,axis=0))
        Ts_A = np.max(Ns,axis=0)-Ns[0,:] + τs[0] 
        Ts_B = np.max(Ns,axis=0)-Ns[-1,:] + τs[-1] 
        f_memory_bell = np.mean(np.exp(-(np.sum(Ts,axis=0)+Ts_A+Ts_B)/τ_coh))

    raw_rate = 1/T_tot
    mu_e2e = mu_link**(2*N_links-1)
    # secret key rate calculations
    f_e2e_qkd = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_qkd
    ex = (1 - mu_e2e)/2
    ez = (1 + mu_e2e)/2 - mu_e2e * f_e2e_qkd
    skr = raw_rate * (1-h([ex])-h([ez]))
    #  fidelity of e2e Bell pairs
    f_e2e_bell = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_bell
    F_e2e = mu_e2e * f_e2e_bell + (1-mu_e2e)/4
    return raw_rate, skr, F_e2e


def T_sequential_no_cutoff(τ_coh, mu_link, F_link,links):
    """ Calculate performance metrics for asynchronous sequential scheme using analytical formulas
    inputs:
        τ_coh: coherence time of quantum memories
        mu_link: parameter in 2qubit depolarizing channel describing noisy link-level entanglement and
        entanglement swapping error
        F_link: fidelity of link level entanglement (i.e.,quality of locally generated Bell pairs)
        links: list of segment (link) lengths in km
    outputs:
        Raw_rate: 1/ expected value of total time for e2e entanglement delivery
        *** application specific quantities:
        skr: secret key rate for qkd (does not include idle times of end memories)
        F_e2e: e2e entanglement fidelity for entanglement distrubtion (does include idle times of end memories)
    """
    if type(links) != np.ndarray:
        links = np.array(links)
    τs = links/c
    T_tot = 2* np.sum( τs / (p_link*Trans(links)) )

    raw_rate = 1/T_tot
    N_links = len(links) # number of links, i.e. no. of repeaters + 1
    mu_e2e = mu_link**(2*N_links-1)
    # secret key rate calculations
    f_memory_qkd = np.prod( p_link*Trans(links[1:])*np.exp(-4*τs[1:]/τ_coh)/(1- (1-p_link*Trans(links[1:]))*np.exp(-2*τs[1:]/τ_coh) )  )
    f_e2e_qkd = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_qkd
    ex = (1 - mu_e2e)/2
    ez = (1 + mu_e2e)/2 - mu_e2e * f_e2e_qkd
    skr = raw_rate * (1-h([ex])-h([ez]))
    #  fidelity of e2e Bell pairs
    Le2e = np.sum(links)
    τe2e = Le2e/c
    f_memory_bell = np.exp(-3*τe2e/τ_coh) *np.prod(p_link*Trans(links[1:])*np.exp(-4*τs[1:]/τ_coh)/(1- (1-p_link*Trans(links[1:]))*np.exp(-4*τs[1:]/τ_coh) ) )
    f_e2e_bell = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_bell
    F_e2e = mu_e2e * f_e2e_bell + (1-mu_e2e)/4

    return raw_rate, skr, F_e2e

def T_parallel_cutoff_slow(τ_cut, τ_coh, mu_link, F_link,links, cct= True, Nmax=100000):
    """ Runs Monte-Carlo iterations to calculate performance metrics for asynchronous parallel scheme
    inputs:
        τ_coh: coherence time of quantum memories
        mu_link: parameter in 2qubit depolarizing channel describing noisy link-level entanglement and
        entanglement swapping error
        F_link: fidelity of link level entanglement (i.e.,quality of locally generated Bell pairs)
        links: list of segment (link) lengths in km
        cct: bool flag to turn on/off classical comm.
        Nmax: No. of MC iterations (if you see oscillations in secret key curves, you need to try larger numbers)

    outputs:
        Raw_rate: 1/ expected value of total time for e2e entanglement delivery
        *** application specific quantities:
        skr: secret key rate for qkd (does not include idle times of end memories)
        F_e2e: e2e entanglement fidelity for entanglement distrubtion (does include idle times of end memories)
    """
    Nmax = int(Nmax) # to make sure Nmax is an integer
    N_links = len(links) # number of links, i.e. no. of repeaters + 1
    if type(links) != np.ndarray:
        links = np.array(links)
    τs = links/c
    Ns = np.zeros((N_links,Nmax))
    for i in range(N_links):
        # Ns[i,:] = F_geo(np.random.rand(Nmax),p_link*Trans(links[i]))
        p_i = p_link*Trans(links[i])
        Ns[i,:] = F_geo(np.random.uniform(low=p_i, high=1, size=(Nmax,)),p_i)

    Tswap = np.zeros((N_links-1,Nmax))
    Tswap_cc = np.zeros((N_links-1,Nmax))
    tL = np.zeros((N_links-1,Nmax))
    tR = np.zeros((N_links-1,Nmax))
    for i in range(N_links-1):
        Tswap[i,:] = np.max(np.stack((2*Ns[i+1,:]*τs[i+1],(2*Ns[i,:]-1)*τs[i]),axis=0), axis=0)
        tL[i,:] = Tswap[i,:]- (2*Ns[i,:]-1)*τs[i]
        tR[i,:] = Tswap[i,:]- 2*(Ns[i+1,:]-1)*τs[i+1]
        Tswap_cc[i,:] = Tswap[i,:] + np.sum(τs[:(i+1)])

    indsL = np.argwhere(np.sum( tL <= τ_cut , axis = 0)== N_links-1 )[:,0]
    indsR = np.argwhere(np.sum( tR <= τ_cut , axis = 0)== N_links-1 )[:,0]
    succ_inds = np.intersect1d(indsL,indsR)
    T_succ = np.sum(np.max(Tswap_cc[:,succ_inds],axis=0))           
    mean_succ =  T_succ/len(succ_inds)
    # return mean_succ
    fail_inds = list(set(list(range(Nmax)))-set(succ_inds))
    T_fail = 0
    for i_r in fail_inds:
        indsL = np.argwhere( tL[:,i_r] > τ_cut )[:,0]
        indsR = np.argwhere( tR[:,i_r] > τ_cut )[:,0]
        tsL = []
        for i_L in indsL:
            tsL.append((2*Ns[i_L,i_r]-1)*τs[i_L] + np.sum(τs[:(i_L+1)]))
        tsR = []
        for i_R in indsR:
            tsR.append(2*(Ns[i_R+1,i_r]-1)*τs[i_R+1]+ np.sum(τs[:(i_R+1)]))
        
        T_fail += min(tsL + tsR)+ τ_cut
        
    raw_rate = len(succ_inds)/(T_succ+T_fail)
    # mean_fail= 0
    # if len(fail_inds)> 0:
    #     mean_fail=  T_fail/len(fail_inds)
    Ts_A = np.max(Tswap_cc[:,succ_inds],axis=0)- 2*(Ns[0,succ_inds]-1)* τs[0] 
    Ts_B = np.max(Tswap_cc[:,succ_inds],axis=0)- (2*Ns[-1,succ_inds]-1)* τs[-1] 
    f_memory_bell = np.mean(np.exp(-(np.sum(tL[:,succ_inds]+tR[:,succ_inds],axis=0)+Ts_A+Ts_B)/τ_coh))

    f_memory_qkd = np.mean(np.exp(-np.sum(tL[:,succ_inds]+tR[:,succ_inds],axis=0)/τ_coh))
    mu_e2e = mu_link**(2*N_links-1)
    # secret key rate calculations
    f_e2e_qkd = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_qkd
    ex = (1 - mu_e2e)/2
    ez = (1 + mu_e2e)/2 - mu_e2e * f_e2e_qkd
    skr = raw_rate * (1-h([ex])-h([ez]))
    #  fidelity of e2e Bell pairs
    f_e2e_bell = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_bell
    F_e2e = mu_e2e * f_e2e_bell + (1-mu_e2e)/4
    return raw_rate, skr, F_e2e



def T_parallel_cutoff(τ_cut, τ_coh, mu_link, F_link,links, cct= True, Nmax=100000):
    """ Runs Monte-Carlo iterations to calculate performance metrics for asynchronous parallel scheme
    inputs:
        τ_coh: coherence time of quantum memories
        mu_link: parameter in 2qubit depolarizing channel describing noisy link-level entanglement and
        entanglement swapping error
        F_link: fidelity of link level entanglement (i.e.,quality of locally generated Bell pairs)
        links: list of segment (link) lengths in km
        cct: bool flag to turn on/off classical comm.
        Nmax: No. of MC iterations (if you see oscillations in secret key curves, you need to try larger numbers)

    outputs:
        Raw_rate: 1/ expected value of total time for e2e entanglement delivery
        *** application specific quantities:
        skr: secret key rate for qkd (does not include idle times of end memories)
        F_e2e: e2e entanglement fidelity for entanglement distrubtion (does include idle times of end memories)
    """
    Nmax = int(Nmax) # to make sure Nmax is an integer
    N_links = len(links) # number of links, i.e. no. of repeaters + 1
    if type(links) != np.ndarray:
        links = np.array(links)
    τs = links/c
    Ns = np.zeros((N_links,Nmax))
    for i in range(N_links):
        # Ns[i,:] = F_geo(np.random.rand(Nmax),p_link*Trans(links[i]))
        p_i = p_link*Trans(links[i])
        Ns[i,:] = F_geo(np.random.uniform(low=p_i, high=1, size=(Nmax,)),p_i)

    Tswap = np.zeros((N_links-1,Nmax))
    Tswap_cc = np.zeros((N_links-1,Nmax))
    tL = np.zeros((N_links-1,Nmax))
    tR = np.zeros((N_links-1,Nmax))
    for i in range(N_links-1):
        Tswap[i,:] = np.max(np.stack((2*Ns[i+1,:]*τs[i+1],(2*Ns[i,:]-1)*τs[i]),axis=0), axis=0)
        tL[i,:] = Tswap[i,:]- (2*Ns[i,:]-1)*τs[i]
        tR[i,:] = Tswap[i,:]- 2*(Ns[i+1,:]-1)*τs[i+1]
        Tswap_cc[i,:] = Tswap[i,:] + np.sum(τs[:(i+1)])

    indsL = np.argwhere(np.sum( tL <= τ_cut , axis = 0)== N_links-1 )[:,0]
    indsR = np.argwhere(np.sum( tR <= τ_cut , axis = 0)== N_links-1 )[:,0]
    succ_inds = np.intersect1d(indsL,indsR)
    T_succ = np.sum(np.max(Tswap_cc[:,succ_inds],axis=0))           
    # mean_succ =  T_succ/len(succ_inds)
    # return mean_succ
    fail_inds = list(set(list(range(Nmax)))-set(succ_inds))

    comm_time = np.array([np.sum(τs[:(i_r+1)]) for i_r in range(N_links-1)]).reshape((N_links-1, 1))
    tsL = (tL[:,fail_inds] > τ_cut)* ( ( (2*Ns[:-1,fail_inds]-1) * τs[:-1].reshape((N_links-1, 1)) ) + comm_time )
    tsR = (tR[:,fail_inds] > τ_cut)* ( (2*(Ns[1:,fail_inds]-1) * τs[1:].reshape((N_links-1, 1)) ) + comm_time )
    ts = np.concatenate( (tsL,tsR), axis=0)
    ts[ts==0] = inf
    T_fail = np.sum( np.min( ts , axis=0 ) +  τ_cut)
    raw_rate = len(succ_inds)/(T_succ+T_fail)
    # mean_fail= 0
    # if len(fail_inds)> 0:
    #     mean_fail=  T_fail/len(fail_inds)

    Ts_A = np.max(Tswap_cc[:,succ_inds],axis=0)- 2*(Ns[0,succ_inds]-1)* τs[0] 
    Ts_B = np.max(Tswap_cc[:,succ_inds],axis=0)- (2*Ns[-1,succ_inds]-1)* τs[-1] 
    f_memory_bell = np.mean(np.exp(-(np.sum(tL[:,succ_inds]+tR[:,succ_inds],axis=0)+Ts_A+Ts_B)/τ_coh))

    f_memory_qkd = np.mean(np.exp(-np.sum(tL[:,succ_inds]+tR[:,succ_inds],axis=0)/τ_coh))
    mu_e2e = mu_link**(2*N_links-1)
    # secret key rate calculations
    f_e2e_qkd = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_qkd
    ex = (1 - mu_e2e)/2
    ez = (1 + mu_e2e)/2 - mu_e2e * f_e2e_qkd
    skr = raw_rate * (1-h([ex])-h([ez]))
    #  fidelity of e2e Bell pairs
    f_e2e_bell = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_bell
    F_e2e = mu_e2e * f_e2e_bell + (1-mu_e2e)/4
    return raw_rate, skr, F_e2e


def T_sequential_cutoff(τ_cut,τ_coh, mu_link, F_link,links):
    """ Calculate performance metrics for asynchronous sequential scheme using analytical formulas
    inputs:
        τ_cut: cut-off time
        τ_coh: coherence time of quantum memories
        mu_link: parameter in 2qubit depolarizing channel describing noisy link-level entanglement and
        entanglement swapping error
        F_link: fidelity of link level entanglement (i.e.,quality of locally generated Bell pairs)
        links: list of segment (link) lengths in km
    outputs:
        Raw_rate: 1/ expected value of total time for e2e entanglement delivery
        *** application specific quantities:
        skr: secret key rate for qkd (does not include idle times of end memories)
        F_e2e: e2e entanglement fidelity for entanglement distrubtion (does include idle times of end memories)
    """
    if type(links) != np.ndarray:
        links = np.array(links)
    τs = links/c
    # implementing the recursion relation :
    ### Tn = Tn-1 / Pn + ( (1/Pn -1) τ_cut + Nm(ms[n],ps[n])*2*τs[n]/Pn)
    p1 = p_link*Trans(links[0])
    T_tot = 2*τs[0]/p1
    for i_l in np.arange(1,len(links)):
        L = links[i_l]
        m_n = int(τ_cut/(2*τs[i_l]))
        p_n = p_link*Trans(L)
        Nm = lambda x: (1-(1+m_n*x)*(1-x)**m_n)/x
        Pm = 1- (1-p_n)**m_n
        T_tot = T_tot / Pm +  (1/Pm -1)*τ_cut + Nm(p_n)*2*τs[i_l]/Pm

    raw_rate = 1/T_tot
    N_links = len(links) # number of links, i.e. no. of repeaters + 1
    mu_e2e = mu_link**(2*N_links-1)
    # secret key rate calculations
    m_arr = np.floor(τ_cut/(2*τs))
    Pm_arr = 1- (1-p_link*Trans(links))**m_arr
    f_memory_qkd = np.prod( p_link*Trans(links[1:])/Pm_arr[1:] * np.exp(-4*τs[1:]/τ_coh) * (1- (1-p_link*Trans(links[1:]))**m_arr[1:] *np.exp(-2*m_arr[1:]*τs[1:]/τ_coh) )  /(1- (1-p_link*Trans(links[1:]))*np.exp(-2*τs[1:]/τ_coh) )  )
    f_e2e_qkd = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_qkd
    ex = (1 - mu_e2e)/2
    ez = (1 + mu_e2e)/2 - mu_e2e * f_e2e_qkd
    skr = raw_rate * (1-h([ex])-h([ez]))
    #  fidelity of e2e Bell pairs
    Le2e = np.sum(links)
    τe2e = Le2e/c
    f_memory_bell = np.exp(-3*τe2e/τ_coh) *np.prod( p_link*Trans(links[1:])/Pm_arr[1:] * np.exp(-4*τs[1:]/τ_coh) * (1- (1-p_link*Trans(links[1:]))**m_arr[1:] *np.exp(-4*m_arr[1:]*τs[1:]/τ_coh) )  /(1- (1-p_link*Trans(links[1:]))*np.exp(-4*τs[1:]/τ_coh) )  )
    f_e2e_bell = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_bell
    F_e2e = mu_e2e * f_e2e_bell + (1-mu_e2e)/4

    # skr, F_e2e = 0, 0
    return raw_rate, skr, F_e2e

def one_repeater_parallel_cutoff(L1,L2,τ_cut,τ_coh, mu_link, F_link, Nmax=100000):
    """ calculates the performance of parallel scheme with one repeater
    inputs:
        L1,L2: elementary link lengths *** works for arbitrary L1,L2 ****
        τ_cut: cut-off time
        τ_coh: coherence time of quantum memories
        mu_link: parameter in 2qubit depolarizing channel describing noisy link-level entanglement and
        entanglement swapping error
        F_link: fidelity of link level entanglement (i.e.,quality of locally generated Bell pairs)
        Nmax: ensemble size for averaging
    outputs:
        Raw_rate: 1/ expected value of total time for e2e entanglement delivery
        *** application specific quantities:
        skr: secret key rate for qkd (does not include idle times of end memories)
        F_e2e: e2e entanglement fidelity for entanglement distrubtion (does include idle times of end memories)
    """
    Nmax = int(Nmax) # to make sure Nmax is an integer
    τ1 = L1/c
    τ2 = L2/c
    p1 = p_link*Trans(L1)
    p2 = p_link*Trans(L2)
    N1 = F_geo(np.random.uniform(low=p1, high=1, size=(Nmax,)),p1)
    N2 = F_geo(np.random.uniform(low=p2, high=1, size=(Nmax,)),p2)
    # N1 = F_geo(np.random.rand(Nmax),p1)
    # N2 = F_geo(np.random.rand(Nmax),p2)
    Ts = np.max(np.array([(2*N1-1)*τ1,2*N2*τ2]),axis=0) # swap moment
    t1L = Ts-(2*N1-1)*τ1 # elapsed time of leftside repeater's memory 
    t1R = Ts-2*(N2-1)*τ2 # elapsed time of rightside repeater's memory 
    indsL = np.argwhere( t1L <= τ_cut )[:,0]
    indsR = np.argwhere( t1R <= τ_cut )[:,0]
    succ_inds = np.intersect1d(indsL,indsR)
    T_succ = np.sum(Ts[succ_inds]+ τ1) 
    fail_inds = list(set(list(range(Nmax)))-set(succ_inds))
    L_runs = np.argwhere( (2*N1[fail_inds]-1)*τ1 <= 2*N2[fail_inds]*τ2 )[:,0]
    R_runs = np.argwhere( (2*N1[fail_inds]-1)*τ1 > 2*N2[fail_inds]*τ2 )[:,0]
    left_mem_on_idx = list(np.array(fail_inds)[L_runs])
    right_mem_on_idx = list(np.array(fail_inds)[R_runs])
    T_elapsed_left = 2*N1[left_mem_on_idx]*τ1 + τ_cut
    T_elapsed_right = 2*(N2[right_mem_on_idx]-1)*τ2 +τ1+ τ_cut
    T_fail = np.sum(T_elapsed_left)+np.sum(T_elapsed_right)

    raw_rate = len(succ_inds)/(T_succ+T_fail)
    t_idle_qkd = t1L + t1R
    f_memory_qkd = np.mean(np.exp(- t_idle_qkd[succ_inds]/τ_coh) )
    N_links = 2
    mu_e2e = mu_link**(2*N_links-1)
    # secret key rate calculations
    f_e2e_qkd = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_qkd
    ex = (1 - mu_e2e)/2
    ez = (1 + mu_e2e)/2 - mu_e2e * f_e2e_qkd
    skr = raw_rate * (1-h([ex])-h([ez]))
    #  fidelity of e2e Bell pairs
    tA = Ts+ τ1 -2*(N1-1)*τ1 # elapsed time of A's memory 
    tB = Ts+ τ1 -(2*N2-1)*τ2 # elapsed time of B's memory 
    t_idle_bell = tA + tB + t1L + t1R
    mean_sender = tA[succ_inds].mean()
    mean_receiver = tB[succ_inds].mean()
    f_memory_bell = np.mean(np.exp(- t_idle_bell[succ_inds]/τ_coh) )
    f_e2e_bell = 0.5 + 0.5 * (2*F_link-1)**N_links *f_memory_bell
    F_e2e = mu_e2e * f_e2e_bell + (1-mu_e2e)/4
    
    return raw_rate, skr, f_e2e_bell
