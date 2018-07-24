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
    elif r > s:
        return - pairing_v_as(p,q,s,r,g)    
    elif math.floor(p/2) == math.floor(q/2) and math.floor(r/2) == math.floor(s/2) and p != q and r != s:
        return -g/2
    else:
        return 0

g = 1#.003
no_of_states = 12
no_of_prt = 6
E_C = 100


f_h = np.zeros(no_of_prt)
f_p = np.zeros(no_of_states-no_of_prt)
for i in range(len(f_h)):
    f_h[i] = sp_f(i,g,fermi_level_number(no_of_prt))
for a in range(len(f_p)):
    f_p[a] = sp_f(a + no_of_prt,g,fermi_level_number(no_of_prt))

v_pp_hh = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
v_hh_pp = np.zeros(shape=(no_of_prt, no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt)) #can be removed when changing v_hh_pp to v_pp_hh in einsum (use correct indices)
v_hh_hh = np.zeros(shape=(no_of_prt, no_of_prt, no_of_prt, no_of_prt))
v_pp_pp = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt))
for matrix in [v_pp_hh, v_hh_hh, v_pp_pp, v_hh_pp]: #v_hh_pp can be removed when changing v_hh_pp to v_pp_hh in einsum (use correct indices)
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
t_2 = v_pp_hh/f_sign_sum

h_bar = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))


def give_next_t_2(t_2_r):        
    intermediate_big_1 = np.einsum('klcd,dblj->kbcj',v_hh_pp,t_2_r)   
    intermediate_small_1 = np.einsum('klcd,cdik->li',v_hh_pp,t_2_r)    
    intermediate_big_2 = np.einsum('klcd,abkl->abcd',v_hh_pp,t_2_r)   
    intermediate_small_2 = np.einsum('klcd,ackl->ad',v_hh_pp,t_2_r)    
    
    h_bar = v_pp_hh \
        + np.einsum('b,abij->abij',f_p,t_2_r) - np.einsum('a,baij->abij',f_p,t_2_r) \
        - np.einsum('j,abij->abij',f_h,t_2_r) + np.einsum('i,abji->abij',f_h,t_2_r) \
        + 1/2 * np.einsum('abcd,cdij->abij',v_pp_pp,t_2_r) \
        + 1/2 * np.einsum('klij,abkl->abij',v_hh_hh,t_2_r) \
        + 1/2 * ( np.einsum('kbcj,acik->abij',intermediate_big_1,t_2_r) - np.einsum('kacj,bcik->abij',intermediate_big_1,t_2_r) - np.einsum('kbci,acjk->abij',intermediate_big_1,t_2_r) + np.einsum('kaci,bcjk->abij',intermediate_big_1,t_2_r) ) \
        + 1/2 * ( np.einsum('li,ablj->abij',intermediate_small_1,t_2_r) - np.einsum('lj,abli->abij',intermediate_small_1,t_2_r) ) \
        + 1/2 * ( np.einsum('ad,dbij->abij',intermediate_small_2,t_2_r) - np.einsum('bd,daij->abij',intermediate_small_2,t_2_r) ) \
        + 1/4 * np.einsum('abcd,cdij->abij',intermediate_big_2,t_2_r) 
        
    delta_t_2_r = h_bar / f_sign_sum #element-wise division
    return (t_2_r + delta_t_2_r)

def recursion(accuracy,maxit,t_2_r=t_2,E_C_r=E_C,it_no=0):
    E_C_old = E_C_r
    t_2_r = give_next_t_2(t_2_r)
    E_C_r = 1/4 * np.einsum('ijab,abij->',v_hh_pp,t_2_r)
    if abs(E_C_r - E_C_old) < accuracy:
        print("After ", it_no , "recursion steps the following result was reached:" , E_C_r)
        print("The residuum is" , E_C_r - E_C_old)
    elif it_no > maxit:
        print("Number of iterations exceeds set value maxit = ", maxit)
        print("Unconverged result is:", E_C_r)
    else: 
        it_no += 1
        recursion(accuracy,maxit,t_2_r=t_2_r,E_C_r=E_C_r,it_no=it_no)
        
#give_next_t_2(t_2)        
    
recursion(0.00001,25)   


 
