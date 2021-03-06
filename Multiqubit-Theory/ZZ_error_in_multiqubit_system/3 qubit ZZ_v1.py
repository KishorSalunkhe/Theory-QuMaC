# Computing the ZZ errors of three transversely coupled qubits
#given frequencies, bare coupling and anharmonicities
"""
Created on Tue May 20 22:15:43 2020

@author: Sumeru
"""
#=========================================================================================
from qutip import *
import numpy as np
no_of_qubits = 3
#=========================================================================================
def binomial(n, k): #calculates binomial coefficients
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0
#=========================================================================================
def threequbit_h(dim, w1, w2, w3, del1, del2, del3, g12, g13, g23):
    #this function defines the two qubit hamiltonian with delta being anharmonicities
    #the exchange coupling is given by g12
    a1 = tensor(destroy(dim), qeye(dim), qeye(dim))
    a1dag = tensor(create(dim), qeye(dim), qeye(dim))
    a2 = tensor(qeye(dim), destroy(dim), qeye(dim))
    a2dag = tensor(qeye(dim), create(dim), qeye(dim))
    a3 = tensor(qeye(dim), qeye(dim), destroy(dim))
    a3dag = tensor(qeye(dim), qeye(dim), create(dim))
    H1 = w1*(a1dag*a1) + (del1/2)*((a1dag*a1)*(a1dag*a1)-(a1dag*a1))
    H2 =  w2*(a2dag*a2) + (del2/2)*((a2dag*a2)*(a2dag*a2)-(a2dag*a2))
    H3 =  w3*(a3dag*a3) + (del3/2)*((a3dag*a3)*(a3dag*a3)-(a3dag*a3))
    Hc = g12*(a1dag*a2 + a1*a2dag) + g13*(a1dag*a3 + a1*a3dag) + g23*(a2dag*a3 + a2*a3dag)
    Htot = H1 + H2 + H3 + Hc
    return Htot
#=========================================================================================
# The ordering of states in qutip is in number_base(n) format as in 
# {000, 001, 002, ..., 00d, 010, 011, ..., 01d, ..., 0d0, ..., 0dd, ..., ddd}
# We shall first truncate the Hamiltonian to limit the maximum number of excitations
# Next we rearrange the Hamiltonian in the binomial pyramid format given below:
# {0(00), 0(01), 0(10), 1(00), 0(02), 0(11), 0(20), 1(01), 1(10), 2(00), ... ,
# 0(0 d), ..., d(00)}
#=========================================================================================    
def red_h(dim, w1, w2, w3, del1, del2, del3, g12, g13, g23):
    d_hs = binomial(dim + no_of_qubits - 1, no_of_qubits ) #dimension of Hilbert Space
    y = np.zeros(d_hs, dtype =int)
    l = 0
    for i in range(0, dim):
        for j in range(0, i+1):
            for k in range(0, i-j+1):
                y[l] = (dim**2)*j + dim*k + (i - j - k)
                l = l+1
    empty2d = np.zeros([d_hs, d_hs], dtype = complex)
    base_ham = threequbit_h(dim, w1, w2, w3, del1, del2, del3, g12, g13, g23)
    for i1 in range(0, d_hs):
        for j1 in range(0, d_hs):
            empty2d[i1, j1] = base_ham[y[i1], y[j1]]
    return empty2d
#=========================================================================================
"Next task is to diagonalize the Hamiltonian"
"Using Q.eigenenergies function in QuTiP or np.linalg does not arranges the values properly"
"we would like to keep the same order as the original basis"
#=========================================================================================
def zz_khz(dim, w1, w2, w3, del1, del2, del3, g12, g13, g23):
    val, vec = np.linalg.eig(red_h(dim, w1, w2, w3, del1, del2, del3, g12, g13, g23))
    vec_arr = np.argmax(vec, axis = 0) #eigenvectors along columns
    val_sor_ind = np.argsort(vec_arr)
    val_incr = val[val_sor_ind]
    zz_split_tot_12 = np.real(val_incr[8] - val_incr[2] - val_incr[3] + val_incr[0])
    zz_split_tot_13 = np.real(val_incr[7] - val_incr[3] - val_incr[1] + val_incr[0])
    zz_split_tot_23 = np.real(val_incr[5] - val_incr[1] - val_incr[2] + val_incr[0])
    zz_err= np.array([zz_split_tot_12/2, zz_split_tot_13/2, zz_split_tot_23/2])
    print('zz_error in kHz, zz_12, zz_13, zz_23')
    return zz_err*1e6
#=========================================================================================

   
    