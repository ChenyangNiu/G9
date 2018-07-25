#this file can do CCD for the pairing problem
import numpy as np
import math



def fermi_level_number(no_of_prt):
    return no_of_prt - 1


def sp_f(p,g,fermi_level):
    """gives back diagonal Fock matrix element""" #Equation 1.4 of CCM-minted.pdf
    d = 1
    if p > fermi_level: #unoccupied states (particles)
        return math.floor( p/2 )*d
    else: #occupied states (holes)
        return math.floor( p/2 )*d - g/2  


def pairing_v_as(p,q,r,s,g):
    """gives back antisymmetric matrix element from pairing interaction"""
    if p > q: #Antisymmetry
        return - pairing_v_as(q,p,r,s,g)
    elif r > s: #Antisymmetry
        return - pairing_v_as(p,q,s,r,g)    
    elif math.floor(p/2) == math.floor(q/2) and math.floor(r/2) == math.floor(s/2) and p != q and r != s: #Matrix element = -g/2 only for pairs
        return -g/2
    else:
        return 0


def matrix_setup(no_of_prt,no_of_states,g):
    """sets up all relevant matrices and quantities"""
    f_h = np.zeros(no_of_prt)   
    f_p = np.zeros(no_of_states-no_of_prt)
    for i in range(len(f_h)):
        f_h[i] = sp_f(i,g,fermi_level_number(no_of_prt))
    for a in range(len(f_p)):
        f_p[a] = sp_f(a + no_of_prt,g,fermi_level_number(no_of_prt))

    v_hh_pp = np.zeros(shape=(no_of_prt, no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt)) 
    v_hh_hh = np.zeros(shape=(no_of_prt, no_of_prt, no_of_prt, no_of_prt))
    v_pp_pp = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_states-no_of_prt))
    for matrix in [v_hh_hh, v_pp_pp, v_hh_pp]: 
        for p in range(np.size(matrix,0)):
            for q in range(np.size(matrix,1)):
                for r in range(np.size(matrix,2)):
                    for s in range(np.size(matrix,3)):
                        matrix[p,q,r,s] = pairing_v_as(p,q,r,s,g) #filling matrix with matrix elements

    f_sign_sum = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
    for a in range(np.size(f_sign_sum,0)):
        for b in range(np.size(f_sign_sum,1)):
            for i in range(np.size(f_sign_sum,2)):
                for j in range(np.size(f_sign_sum,3)):
                    f_sign_sum[a,b,i,j] = - f_p[a] - f_p[b] + f_h[i] + f_h[j] #this combination of f_p and f_h is needed in some cases

    t_2 = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
    t_2 = v_hh_pp.T/f_sign_sum #initial guess for t_2
    
    return f_h, f_p, v_hh_pp, v_hh_hh, v_pp_pp, f_sign_sum, t_2


def give_next_t_2(t_2, f_h, f_p, v_hh_pp, v_hh_hh, v_pp_pp, f_sign_sum):  
    """gives back the resulting t_2 from one iteration""" #Eq. 1.36 from CCM-minted.pdf
    intermediate_big_1 = np.einsum('klcd,dblj->kbcj',v_hh_pp,t_2) #Intermediates for more efficient calculation of h_bar
    intermediate_small_1 = np.einsum('klcd,cdik->li',v_hh_pp,t_2)    
    intermediate_big_2 = np.einsum('klcd,abkl->abcd',v_hh_pp,t_2)   
    intermediate_small_2 = np.einsum('klcd,ackl->ad',v_hh_pp,t_2)    
    
    #Next lines contain the unnumbered equation on top of p. 25 in CCM-minted.pdf 
    h_bar = v_hh_pp.T \
        + np.einsum('b,abij->abij',f_p,t_2) - np.einsum('a,baij->abij',f_p,t_2) \
        - np.einsum('j,abij->abij',f_h,t_2) + np.einsum('i,abji->abij',f_h,t_2) \
        + 1/2 * np.einsum('abcd,cdij->abij',v_pp_pp,t_2) \
        + 1/2 * np.einsum('klij,abkl->abij',v_hh_hh,t_2) \
        + 1/2 * ( np.einsum('kbcj,acik->abij',intermediate_big_1,t_2) - np.einsum('kacj,bcik->abij',intermediate_big_1,t_2) - np.einsum('kbci,acjk->abij',intermediate_big_1,t_2) + np.einsum('kaci,bcjk->abij',intermediate_big_1,t_2) ) \
        + 1/2 * ( np.einsum('li,ablj->abij',intermediate_small_1,t_2) - np.einsum('lj,abli->abij',intermediate_small_1,t_2) ) \
        + 1/2 * ( np.einsum('ad,dbij->abij',intermediate_small_2,t_2) - np.einsum('bd,daij->abij',intermediate_small_2,t_2) ) \
        + 1/4 * np.einsum('abcd,cdij->abij',intermediate_big_2,t_2) 
    
    delta_t_2 = h_bar / f_sign_sum 
    return (t_2 + delta_t_2) #new value of t_2 is old value plus correction delta_t_2


def recursion(accuracy,maxit,t_2,E_C_r,it_no, f_h, f_p, v_hh_pp, v_hh_hh, v_pp_pp, f_sign_sum):
    """calls give_next_t_2 until the resulting energy doesn't change more than specified with accuracy"""
    E_C_old = E_C_r #set E_C_old to last E_C value
    t_2 = give_next_t_2(t_2, f_h, f_p, v_hh_pp, v_hh_hh, v_pp_pp, f_sign_sum) #calculate new iteration for t_2
    E_C_r = 1/4 * np.einsum('ijab,abij->',v_hh_pp,t_2) #calculate new value for E_C
    if abs(E_C_r - E_C_old) < accuracy: #stop recursion when desired accuracy is reached
        print("After ", it_no , "recursion steps the following result was reached:" , E_C_r)
        print("The residuum is" , E_C_r - E_C_old)
        return E_C_r 
    elif it_no > maxit: #stop recursion when set value of maximal iterations is exceeded
        print("Number of iterations exceeds set value maxit = ", maxit)
        print("Unconverged result is:", E_C_r)
    else: #otherwise continue recursion
        it_no += 1
        recursion(accuracy,maxit,t_2,E_C_r,it_no,f_h, f_p, v_hh_pp, v_hh_hh, v_pp_pp, f_sign_sum)
        

def pairing_ccd_main(no_of_prt, no_of_states, g, accuracy, maxit):
    """main function for CCD calculation for the pairing model"""
    f_h, f_p, v_hh_pp, v_hh_hh, v_pp_pp, f_sign_sum, t_2 = matrix_setup(no_of_prt,no_of_states,g)
    E_C = 100    
    res_E_C = recursion(accuracy, maxit, t_2, E_C, 0, f_h, f_p, v_hh_pp, v_hh_hh, v_pp_pp, f_sign_sum)
    print(res_E_C)
    return res_E_C

  


 
