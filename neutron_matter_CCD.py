#this file maybe sometime can do CCD for neutron matter with the Minnesota potential
import numpy as np
import math

#hbarc = 197.32697
A_magic = 14
en_magic = 1 #maximal k**2 for states within A_magic
density = 0.1
grid_spacing = (A_magic / density) ** 1/3
momentum_cutoff = (2*0.465) ** 1/2
n_max = 2

def check_n_max(n_max):
    if 2 * math.pi * n_max / grid_spacing <= momentum_cutoff:
        print("n_max is too small")

def make_sp_states(n_max, A_magic, en_magic):
    grid_dim = 2*(2*n_max+1)**3

    #all_sp_states = np.zeros(shape=(grid_dim,4))
    hole_sp_states = np.zeros(shape=(A_magic,4))
    prt_sp_states = np.zeros(shape=(grid_dim - A_magic,4))    
    m_s = -1
    n_z = -n_max
    n_y = -n_max
    n_x = -n_max
    i = 0
    j = 0
    for p in range(grid_dim): #this fills up hole_sp_states and prt_sp_states with the correct single-particle states
        if math.floor(n_x)**2 + math.floor(n_y)**2 + math.floor(n_z)**2 <= en_magic:
            hole_sp_states[i] = [math.floor(n_x),math.floor(n_y),math.floor(n_z),m_s] 
            i += 1
        else: 
            prt_sp_states[j] = [math.floor(n_x),math.floor(n_y),math.floor(n_z),m_s] 
            j += 1            
        #all_sp_states[p] = [math.floor(n_x),math.floor(n_y),math.floor(n_z),m_s]
        m_s += 2
        if m_s > 1:
            m_s = -1
        n_z += 0.5
        if round(n_z,6) > n_max + 1 - 0.5:
            n_z = -n_max
        n_y += (0.5 / (2*n_max+1))
        if round(n_y,6) > n_max + 1 - (0.5 / (2*n_max+1)):
            n_y = -n_max
        n_x += (0.5 / (2*n_max+1)**2)
        if round(n_x,6) > n_max + 1 - (0.5 / (2*n_max+1)**2):
            n_x = -n_max
    
    #print(hole_sp_states)        
    #print(prt_sp_states)
    return hole_sp_states, prt_sp_states


def make_tp_states(n_max, A_magic, en_magic):
    hole_sp_states, prt_sp_states = make_sp_states(n_max, A_magic, en_magic)
    
    hh_states = np.zeros(shape=( A_magic*(A_magic-1) ,8)) #this is a list of all hh states with center-of-mass and relative momentum
    k = 0
    for i in range(A_magic):
        for j in range(A_magic):
            if j != i:
                for coord in range(3):
                    hh_states[k, coord] = hole_sp_states[i, coord] + hole_sp_states[j, coord] #first 3 entries of hh_states are the center-of-mass momentum
                    hh_states[k, coord+3] = hole_sp_states[i, coord] - hole_sp_states[j, coord] #next 3 entries are 2 * relative momentum
                hh_states[k, 6] = hole_sp_states[i, 3]
                hh_states[k, 7] = hole_sp_states[j, 3]
                k += 1
    
    hh_com_momenta = np.zeros(shape=(0 ,4)) #this a list of center-of-mass momenta of all hh states; the last entry will later contain the number of correponding states
    for n_x in range(-n_max, n_max+1):
        for n_y in range(-n_max, n_max+1):
            for n_z in range(-n_max, n_max+1):
                if abs(n_x) + abs(n_y) + abs(n_z) <= 2*en_magic**(1/2):
                    hh_com_momenta = np.append(hh_com_momenta,[[n_x,n_y,n_z,0]],axis=0)
    
    hh_channel_states = np.zeros(shape=( A_magic*(A_magic-1) ,5)) #this is a list of all hh states with relative momentum, where the states are sorted according to their total momentum (channel)
    l = 0
    for cm in hh_com_momenta:
        old_l = l
        for st in hh_states: 
            if np.array_equal(st[0:3],cm[0:3]):
                hh_channel_states[l] = st[3:8]
                l += 1
                print(l)
        hh_com_momenta[i,3] = l - old_l

    #print(hh_channel_states)             
    print(hh_com_momenta)
    print(len(hh_com_momenta))
                                     
    return hh_com_momenta, hh_channel_states

    
make_tp_states(n_max, A_magic, en_magic)