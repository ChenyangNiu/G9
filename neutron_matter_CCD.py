#this file maybe sometime can do CCD for neutron matter with the Minnesota potential
import numpy as np
import math

hbarc = 197.32697
mass = 939.56541 #neutron mass
A_magic = 14
en_magic = 1 #maximal k**2 (k_Fermi) for states within A_magic
density = 0.1
grid_spacing = (A_magic / density) ** (1/3)
momentum_cutoff = (2*0.465) ** (1/2)
n_max = 2


def check_n_max(n_max):
    if 2 * math.pi * n_max / grid_spacing <= 1.5 * momentum_cutoff:
        print("n_max is too small")
        
   
def grid_dim(n_max):
    return 2*(2*n_max+1)**3
        

def make_sp_states(n_max, A_magic, en_magic):
    """This makes all hole and particle single-body states."""
    hole_sp_states = np.zeros(shape=(A_magic,4), dtype=int)
    prt_sp_states = np.zeros(shape=(grid_dim(n_max) - A_magic,4), dtype=int)    
    m_s = -1
    n_z = -n_max
    n_y = -n_max
    n_x = -n_max
    i = 0
    j = 0
    for p in range(grid_dim(n_max)): #this fills up hole_sp_states and prt_sp_states with the correct single-particle states #yes, I know, it is way too complicated
        if math.floor(n_x)**2 + math.floor(n_y)**2 + math.floor(n_z)**2 <= en_magic:
            hole_sp_states[i] = [math.floor(n_x),math.floor(n_y),math.floor(n_z),m_s] 
            i += 1
        else: 
            prt_sp_states[j] = [math.floor(n_x),math.floor(n_y),math.floor(n_z),m_s] 
            j += 1            
        m_s += 2 #the following can be way easier done with several loops instead of just one
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
    
    return hole_sp_states, prt_sp_states


def make_tp_states(n_max, A_magic, en_magic):
    """This makes all needed hole and particle two-body states."""
    hole_sp_states, prt_sp_states = make_sp_states(n_max, A_magic, en_magic)
    
    #hh:
    hh_states = np.zeros(shape=( A_magic*(A_magic-1) ,8), dtype=int) #this is a list of all hh states with center-of-mass and relative momentum
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
    
    com_momenta = np.zeros(shape=(0 ,7), dtype=int) #this a list of center-of-mass momenta of all hh states; the last entries will later contain the number of correponding states
    for n_x in range(-n_max, n_max+1):
        for n_y in range(-n_max, n_max+1):
            for n_z in range(-n_max, n_max+1):
                if abs(n_x) + abs(n_y) + abs(n_z) <= 2*en_magic**(1/2):
                    com_momenta = np.append(com_momenta,[[n_x,n_y,n_z,0,0,0,0]],axis=0)
    
    hh_channel_states = np.zeros(shape=( A_magic*(A_magic-1) ,5), dtype=int) #this is a list of all hh states with relative momentum, where the states are sorted according to their total momentum (channel)
    l = 0
    for i in range(len(com_momenta)):
        old_l = l
        for st in hh_states: 
            if np.array_equal(st[0:3],com_momenta[i,0:3]):
                hh_channel_states[l] = st[3:8]
                l += 1
        com_momenta[i,3] = round(old_l,0) #position in hh_channel_states where the states corresponding to the center-of-mass momentum in com_momenta[i,0:3] starts
        com_momenta[i,4] = round(l - old_l,0) #number of states in hh_channel_states corresponding to this channel
    
    #pp:
    pp_channel_states = np.zeros(shape=(0 ,5), dtype=int) #this is a list of all pp states with relative momentum, where the states are sorted according to their total momentum (channel)
    l = 0
    for k in range(len(com_momenta)):#for every hole-hole center-of-mass momentum ...
        old_l = l
        for i in range(grid_dim(n_max) - A_magic): #...and every single-particle state...
            for j in range(grid_dim(n_max) - A_magic): #...and every other single-particle state...
                if j != i: 
                    current_pp_com_momentum = prt_sp_states[i, 0:3] + prt_sp_states[j, 0:3] #...check if the single-particle states can form this center-of-mass momentum;...
                    if np.array_equal(current_pp_com_momentum, com_momenta[k,0:3]):#...if that is true, write the corresponding relative momentum in the correct channel   
                        append_array = np.zeros(5)
                        append_array[0:3] = prt_sp_states[i, 0:3] - prt_sp_states[j, 0:3]
                        append_array[3] = prt_sp_states[i, 3]
                        append_array[4] = prt_sp_states[j, 3]    
                        pp_channel_states = np.append(pp_channel_states,[append_array],axis=0)    
                        l += 1       
        com_momenta[k,5] = old_l #position in pp_channel_states where the states corresponding to the center-of-mass momentum in com_momenta[i,0:3] starts
        com_momenta[k,6] = l - old_l #number of states in pp_channel_states corresponding to this channel
    
    #ph:
    ph_com_momenta = np.zeros(shape=(0 ,5), dtype=int) #this a list of center-of-mass momenta of all ph states; the last entries will later contain the number of correponding states    
    ph_states = np.zeros(shape=( A_magic*(grid_dim(n_max) - A_magic) ,8), dtype=int) #this is a list of all ph states with center-of-mass and relative momentum
    k = 0
    for i in range(grid_dim(n_max) - A_magic):
        for j in range(A_magic):
            current_ph_com_momentum = prt_sp_states[i,0:3] + hole_sp_states[j,0:3]
            for coord in range(3):
                ph_states[k, coord] = current_ph_com_momentum[coord] #first 3 entries of ph_states are the center-of-mass momentum
                ph_states[k, coord+3] = prt_sp_states[i, coord] - hole_sp_states[j, coord] #next 3 entries are 2 * relative momentum
            ph_states[k, 6] = prt_sp_states[i, 3]
            ph_states[k, 7] = hole_sp_states[j, 3]
            k += 1
            ph_com_momentum_found = False
            for l in range(len(ph_com_momenta)):
                if np.array_equal(current_ph_com_momentum,ph_com_momenta[l,0:3]):
                    ph_com_momentum_found = True
                    break
            if ph_com_momentum_found == False:
                ph_com_momenta = np.append(ph_com_momenta, [np.concatenate((current_ph_com_momentum,[0,0]))], axis=0)      
                       
    ph_channel_states = np.zeros(shape=( A_magic*(grid_dim(n_max) - A_magic) ,5), dtype=int) #this is a list of all ph states with relative momentum, where the states are sorted according to their total momentum (channel)
    l = 0
    for i in range(len(ph_com_momenta)):
        old_l = l
        for st in ph_states: 
            if np.array_equal(st[0:3],ph_com_momenta[i,0:3]):
                ph_channel_states[l] = st[3:8]
                l += 1
        ph_com_momenta[i,3] = round(old_l,0) #position in ph_channel_states where the states corresponding to the center-of-mass momentum in ph_com_momenta[i,0:3] starts
        ph_com_momenta[i,4] = round(l - old_l,0) #number of states in ph_channel_states corresponding to this channel
        
    print(com_momenta)        
    print(ph_com_momenta)                                   
    return com_momenta, hh_channel_states, pp_channel_states, ph_com_momenta, ph_channel_states


def save_tp_states(n_max, A_magic, en_magic): 
    """This function saves calculated hh and pp states. They can be loaded with load_tp_states, which is much faster than recalculating them again."""
    com_momenta, hh_channel_states, pp_channel_states, ph_com_momenta, ph_channel_states = make_tp_states(n_max, A_magic, en_magic)
    np.save("com_momenta-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy", com_momenta)
    np.save("hh_channel_states-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy", hh_channel_states)
    np.save("pp_channel_states-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy", pp_channel_states)
    np.save("ph_com_momenta-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy", ph_com_momenta)
    np.save("ph_channel_states-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy", ph_channel_states)


def load_tp_states(n_max, A_magic, en_magic):
    com_momenta = np.load("com_momenta-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy")
    hh_channel_states = np.load("hh_channel_states-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy")
    pp_channel_states = np.load("pp_channel_states-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy")
    ph_com_momenta = np.load("ph_com_momenta-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy")
    ph_channel_states = np.load("ph_channel_states-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+".npy")
    return com_momenta, hh_channel_states, pp_channel_states, ph_com_momenta, ph_channel_states
    

def give_minnesota_constants(grid_spacing):
    r_factor = 200 / grid_spacing**3 * (math.pi/1.487)**(3/2)
    s_factor = -91.85 / grid_spacing**3 * (math.pi/0.465)**(3/2)
    r_exp = 1/(4*1.487)
    s_exp = 1/(4*0.465)
    q_square_factor = (math.pi / grid_spacing) ** 2 #1/2 of q cancles with 2 of 2 pi
    f_factor = 2 * hbarc**2 / mass * q_square_factor
    return r_factor, s_factor, r_exp, s_exp, q_square_factor, f_factor


def give_neutron_minnesota_V(bra_rel_numbers,bra_s_1,bra_s_2,ket_rel_numbers,ket_s_1,ket_s_2, r_factor,s_factor, r_exp,s_exp, q_square_factor):
    """Calculate the Minnesota matrix element for two-body states bra and ket."""
    if bra_s_1 == bra_s_2 or ket_s_1 == ket_s_2: #this also accounts for Pauli principle
        return 0
    elif bra_s_1 == ket_s_1: # then automatically bra_s_2 == ket_s_2 holds
        q_mi_square = q_square_factor * np.einsum('i,i->', (bra_rel_numbers - ket_rel_numbers), (bra_rel_numbers - ket_rel_numbers))
        q_pl_square = q_square_factor * np.einsum('i,i->', (bra_rel_numbers + ket_rel_numbers), (bra_rel_numbers + ket_rel_numbers))
        return (r_factor * math.exp(- q_mi_square * r_exp) + s_factor * math.exp(- q_mi_square * s_exp) 
                + r_factor * math.exp(- q_pl_square * r_exp) + s_factor * math.exp(- q_pl_square * s_exp)) #antisymmetrized matrix element
    else: # then automatically bra_s_1 == ket_s_2 and bra_s_2 == ket_s_1 holds
        return - give_neutron_minnesota_V(bra_rel_numbers,bra_s_2,bra_s_1,ket_rel_numbers,ket_s_1,ket_s_2, r_factor,s_factor, r_exp,s_exp, q_square_factor)   


        
def give_sp_f(com_numbers, state_numbers, which_spin_number, r_factor, s_factor, r_exp, s_exp, q_square_factor, f_factor, n_max, A_magic, en_magic, grid_spacing):
    hole_sp_states, prt_sp_states = make_sp_states(n_max, A_magic, en_magic)
    sp_momentum_numbers = (com_numbers + (-1)**which_spin_number * state_numbers[0:3])/2 #this is the resulting single-particle momentum expressed in n's
    kin_energy = f_factor * np.einsum('i,i->', sp_momentum_numbers, sp_momentum_numbers)  
    pot_sum_energy = 0 
    for st in hole_sp_states: 
        pot_sum_energy += give_neutron_minnesota_V(st[0:3] - sp_momentum_numbers, st[3], state_numbers[which_spin_number + 3], st[0:3] - sp_momentum_numbers, st[3], state_numbers[which_spin_number + 3], r_factor,s_factor, r_exp,s_exp, q_square_factor)
    return kin_energy + pot_sum_energy        
    

def matrix_setup(n_max, A_magic, en_magic, grid_spacing):
    """Sets up F, V, and T_2 as matrices.""" 
    com_momenta, hh_channel_states, pp_channel_states, ph_com_momenta, ph_channel_states = load_tp_states(n_max, A_magic, en_magic)#make_tp_states(n_max, A_magic, en_magic)
    r_factor, s_factor, r_exp, s_exp, q_square_factor,f_factor = give_minnesota_constants(grid_spacing)
    no_of_com_momenta = len(com_momenta)
    no_of_ph_com_momenta = len(ph_com_momenta)
    f_h_list = []
    f_p_list = []
    v_pp_hh_list = []
    v_hh_hh_list = []
    v_pp_pp_list = []
    v_ph_ph_list = []
    f_sign_sum_list = []
    t_2_list = []
    
    for i in range(no_of_com_momenta):
        f_h_list.append(np.zeros(shape=(com_momenta[i,4],2)))
        for j in range(com_momenta[i,3],com_momenta[i,3]+com_momenta[i,4]):
            f_h_list[i][j - int(com_momenta[i,3]), 0] = give_sp_f(com_momenta[i,0:3], hh_channel_states[j], 0, r_factor, s_factor, r_exp, s_exp, q_square_factor, f_factor, n_max, A_magic, en_magic, grid_spacing)
            f_h_list[i][j - int(com_momenta[i,3]), 1] = give_sp_f(com_momenta[i,0:3], hh_channel_states[j], 1, r_factor, s_factor, r_exp, s_exp, q_square_factor, f_factor, n_max, A_magic, en_magic, grid_spacing)
            
    for i in range(no_of_com_momenta):
        f_p_list.append(np.zeros(shape=(com_momenta[i,6],2)))
        for j in range(com_momenta[i,5],com_momenta[i,5]+com_momenta[i,6]):
            f_p_list[i][j - int(com_momenta[i,5]), 0] = give_sp_f(com_momenta[i,0:3], pp_channel_states[j], 0, r_factor, s_factor, r_exp, s_exp, q_square_factor, f_factor, n_max, A_magic, en_magic, grid_spacing)
            f_p_list[i][j - int(com_momenta[i,5]), 1] = give_sp_f(com_momenta[i,0:3], pp_channel_states[j], 1, r_factor, s_factor, r_exp, s_exp, q_square_factor, f_factor, n_max, A_magic, en_magic, grid_spacing)
    
    for i in range(no_of_com_momenta):
        v_pp_hh_list.append(np.zeros(shape=(com_momenta[i,6], com_momenta[i,4])))
        for j in range(com_momenta[i,5],com_momenta[i,5]+com_momenta[i,6]):
            for k in range(com_momenta[i,3],com_momenta[i,3]+com_momenta[i,4]):
                v_pp_hh_list[i][j - int(com_momenta[i,5]), k - int(com_momenta[i,3])] = give_neutron_minnesota_V(pp_channel_states[j,0:3],pp_channel_states[j,3],pp_channel_states[j,4], hh_channel_states[k,0:3],hh_channel_states[k,3],hh_channel_states[k,4], r_factor,s_factor, r_exp,s_exp, q_square_factor)
                
        v_hh_hh_list.append(np.zeros(shape=(com_momenta[i,4], com_momenta[i,4])))
        for j in range(com_momenta[i,3],com_momenta[i,3]+com_momenta[i,4]):
            for k in range(com_momenta[i,3],com_momenta[i,3]+com_momenta[i,4]):
                v_hh_hh_list[i][j - int(com_momenta[i,3]), k - int(com_momenta[i,3])] = give_neutron_minnesota_V(hh_channel_states[j,0:3],hh_channel_states[j,3],hh_channel_states[j,4], hh_channel_states[k,0:3],hh_channel_states[k,3],hh_channel_states[k,4], r_factor,s_factor, r_exp,s_exp, q_square_factor)
                
        v_pp_pp_list.append(np.zeros(shape=(com_momenta[i,6], com_momenta[i,6])))
        for j in range(com_momenta[i,5],com_momenta[i,5]+com_momenta[i,6]):
            for k in range(com_momenta[i,5],com_momenta[i,5]+com_momenta[i,6]):
                v_pp_pp_list[i][j - int(com_momenta[i,5]), k - int(com_momenta[i,5])] = give_neutron_minnesota_V(pp_channel_states[j,0:3],pp_channel_states[j,3],pp_channel_states[j,4], pp_channel_states[k,0:3],pp_channel_states[k,3],pp_channel_states[k,4], r_factor,s_factor, r_exp,s_exp, q_square_factor)
            
    #for i in range(no_of_ph_com_momenta):                
        #v_ph_ph_list.append(np.zeros(shape=(ph_com_momenta[i,4], ph_com_momenta[i,4])))
        #for j in range(ph_com_momenta[i,3],ph_com_momenta[i,3]+ph_com_momenta[i,4]):
            #for k in range(ph_com_momenta[i,3],ph_com_momenta[i,3]+ph_com_momenta[i,4]):
                #v_ph_ph_list[i][j - int(ph_com_momenta[i,3]), k - int(ph_com_momenta[i,3])] = give_neutron_minnesota_V(ph_channel_states[j,0:3],ph_channel_states[j,3],ph_channel_states[j,4], ph_channel_states[k,0:3],ph_channel_states[k,3],ph_channel_states[k,4], r_factor,s_factor, r_exp,s_exp, q_square_factor)
        
    for i in range(no_of_com_momenta):
        f_sign_sum_list.append(np.zeros(shape=(com_momenta[i,6], com_momenta[i,4])))
        for j in range(com_momenta[i,5],com_momenta[i,5]+com_momenta[i,6]):
            for k in range(com_momenta[i,3],com_momenta[i,3]+com_momenta[i,4]):
                f_sign_sum_list[i][j - int(com_momenta[i,5]), k - int(com_momenta[i,3])] = f_h_list[i][k - int(com_momenta[i,3]), 0] + f_h_list[i][k - int(com_momenta[i,3]), 1] - f_p_list[i][j - int(com_momenta[i,5]), 0] - f_p_list[i][j - int(com_momenta[i,5]), 1]
    
    for i in range(no_of_com_momenta):
        t_2_list.append([0])
        t_2_list[i] = v_pp_hh_list[i] / f_sign_sum_list[i]
            
    return v_pp_hh_list, v_hh_hh_list, v_pp_pp_list, v_ph_ph_list, f_h_list, f_p_list, f_sign_sum_list, t_2_list   


def save_matrices(n_max, A_magic, en_magic, grid_spacing):
    """This function saves calculated matrices. They can be loaded with load_matrices, which is much faster than recalculating them again."""
    v_pp_hh_list, v_hh_hh_list, v_pp_pp_list, v_ph_ph_list, f_h_list, f_p_list, f_sign_sum_list, t_2_list = matrix_setup(n_max, A_magic, en_magic, grid_spacing)
    np.save("v_pp_hh-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", v_pp_hh_list)
    np.save("v_hh_hh-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", v_hh_hh_list)
    np.save("v_pp_pp-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", v_pp_pp_list)
    #np.save("v_ph_ph-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", v_ph_ph_list)
    np.save("f_h-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", f_h_list)
    np.save("f_p-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", f_p_list)
    np.save("f_sign_sum-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", f_sign_sum_list)
    np.save("t_2-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy", t_2_list)


def load_matrices(n_max, A_magic, en_magic, grid_spacing):
    v_pp_hh_list = np.load("v_pp_hh-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")
    v_hh_hh_list = np.load("v_hh_hh-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")
    v_pp_pp_list = np.load("v_pp_pp-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")   
    #v_ph_ph_list = np.load("v_ph_ph-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")  
    v_ph_ph_list = [] #as it is not needed anywhere
    f_h_list = np.load("f_h-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")
    f_p_list = np.load("f_p-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")
    f_sign_sum_list = np.load("f_sign_sum-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")
    t_2_list = np.load("t_2-n_max_"+str(n_max)+"-A_magic_"+str(A_magic)+"-grid_spacing_"+str(int(grid_spacing))+".npy")
    return v_pp_hh_list, v_hh_hh_list, v_pp_pp_list, v_ph_ph_list, f_h_list, f_p_list, f_sign_sum_list, t_2_list


#

com_momenta, hh_channel_states, pp_channel_states, ph_com_momenta, ph_channel_states = load_tp_states(n_max, A_magic, en_magic)
v_pp_hh_list, v_hh_hh_list, v_pp_pp_list, v_ph_ph_list, f_h_list, f_p_list, f_sign_sum_list, t_2_list = load_matrices(n_max, A_magic, en_magic, grid_spacing)
r_factor, s_factor, r_exp, s_exp, q_square_factor, f_factor = give_minnesota_constants(grid_spacing)

print(com_momenta)
print(len(com_momenta))
print(hh_channel_states[0:2])
print(t_2_list[0])

#save_tp_states(n_max, A_magic, en_magic)
#save_matrices(n_max, A_magic, en_magic, grid_spacing)


