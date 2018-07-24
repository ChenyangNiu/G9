#this file computes the pairing matrix for a given number of pairs and levels and solves the eigenvalue problem
import numpy as np
from numpy import linalg as la
from sympy.utilities.iterables import multiset_permutations as mup
import matplotlib.pyplot as plt
import math
import array

d = 1

def binomial(a,b):
    return math.factorial(a)/(math.factorial(a-b)*math.factorial(b))


def fermi_level_number(no_of_prt):
    return no_of_prt - 1


def sp_f(p,g,fermi_level):
    #Equation 1.4 of CCM-minted.pdf
    if p > fermi_level: #unoccupied states (particles)
        return math.ceil( (p-fermi_level)/2 )*d
    else: #occupied states (holes)
        return math.ceil( (p-fermi_level)/2 )*d - g/2    


def pairing_v_as(p,q,r,s,g):
    if p > q:
        return - pairing_v_as(q,p,r,s,g)
    elif s > r:
        return - pairing_v_as(p,q,s,r,g)    
    elif math.floor(p/2) == math.floor(q/2) and math.floor(r/2) == math.floor(s/2) and p != q and r != s:
        return -g/2
    else:
        return 0


g = .3
no_of_states = 14
no_of_prt = 6


f_h = np.zeros(no_of_prt)
f_p = np.zeros(no_of_states-no_of_prt)
for i in range(len(f_h)):
    f_h[i] = sp_f(i,g,fermi_level_number(no_of_prt))
for a in range(len(f_p)):
    f_p[a] = sp_f(a + no_of_prt,g,fermi_level_number(no_of_prt))

v_pp_hh = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
v_hh_hh = np.zeros(shape=(no_of_prt, no_of_prt, no_of_prt, no_of_prt))
v_pp_pp = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt))
for matrix in [v_pp_hh, v_hh_hh, v_pp_pp]:
    for p in range(np.size(matrix,0)):
        for q in range(np.size(matrix,1)):
            for r in range(np.size(matrix,2)):
                for s in range(np.size(matrix,3)):
                    matrix[p,q,r,s] = pairing_v_as(p,q,r,s,g)

f_sign_sum = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
for a in range(np.size(f_sign_sum,0)):
    for b in range(np.size(f_sign_sum,1)):
        for i in range(np.size(f_sign_sum,2)):
            for j in range(np.size(f_sign_sum,3)):
                f_sign_sum[a,b,i,j] = f_p[a] + f_p[b] - f_h[i] - f_h[j]

t_2 = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
for a in range(np.size(t_2,0)):
    for b in range(np.size(t_2,1)):
        for i in range(np.size(t_2,2)):
            for j in range(np.size(t_2,3)):
                t_2[a,b,i,j] = v_pp_hh[a,b,i,j]/f_sign_sum[a,b,i,j]

h_bar = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
                
#print(f_h)
#print(t_2[2,3,:,:])
#print(v_pp_hh[2,3,:,:])

def give_next_t_2(t_2):    
    h_bar = v_pp_hh \
        + np.einsum('b,abij->abij',f_p,t_2) - np.einsum('a,baij->abij',f_p,t_2) \
        - np.einsum('j,abij->abij',f_p,t_2) + np.einsum('i,abji->abij',f_p,t_2)
    delta_t_2 = h_bar / f_sign_sum #element-wise division
    return t_2 + delta_t_2

print(h_bar[1,0])
print(f_p[0]*t_2[1,0,:,:])




#print(V[2,3,no_of_prt:,no_of_prt:])
#print(t[0,1])

#fl = fermi_level_number(no_of_prt)    
#testest = np.einsum('cd,cd->',V[2,3,no_of_prt:,no_of_prt:],t[0,1])
#print(testest)
 
