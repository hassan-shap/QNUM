import numpy as np
import matplotlib.pyplot as plt
from math import *


def ent_rate_memory(p_list,W):
    N = len(p_list) # number of links (or number of repeaters -1)
    Q = np.zeros((W+1,N+1)) 
    for k in range(1,N+1):
        p_k = p_list[k-1]
        for i in range(W+1):
            Q[i,k] = comb(W,i)* p_k**(i) * (1-p_k)**(W-i)
            
    P = np.zeros((W+1,N+1)) 
    P[:,1] = Q[:,1]     
    for k in range(2,N+1):
        for i in range(W+1):
            P[i,k] = P[i,k-1]* np.sum(Q[i:,k]) + Q[i,k]* np.sum(P[i+1:,k-1])

    return np.sum(np.arange(W+1)*P[:,N])

    
q = 1.0
p = 0.9
N_list = np.arange(1,11)
W_list = [1,2,3]
Et = np.zeros((len(W_list),len(N_list)))
plt.figure(figsize=(5,3))
for i_W, W in enumerate(W_list):
    for i_N, N in enumerate(N_list):
        p_list = [p]*N
        Et[i_W,i_N] = q**(N-1) * ent_rate_memory(p_list,W)
    plt.plot(N_list, Et[i_W,:],".-",label="W=%d" % W)
plt.legend()
plt.grid()
plt.ylabel("Expected throughput")
plt.xlabel("Number of hops")
plt.title("p=%.1f" % p)
plt.show()