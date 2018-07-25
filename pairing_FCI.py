#this file computes the pairing matrix for a given number of pairs and levels and solves the eigenvalue problem
import numpy as np
from numpy import linalg as la
from sympy.utilities.iterables import multiset_permutations as mup
import matplotlib.pyplot as plt
import math
import array


def binomial(a,b):
    return math.factorial(a)/(math.factorial(a-b)*math.factorial(b))


def make_slater_determinants(no_of_pairs,no_of_levels):
    no_of_slater_determinants = int(binomial(no_of_levels, no_of_pairs)) #number of Slater determinants
    creation_array = np.concatenate((np.ones(no_of_pairs), np.zeros(no_of_levels-no_of_pairs))) #Array with ones 
    #for pairs, other levels are filled with zeros
    slater_determinants = np.zeros(shape=(no_of_slater_determinants,no_of_levels))
    i = 0
    for p in mup(creation_array):
        slater_determinants[i] = p #Slater determinants are set as all distinct permutations of the creationarray
        i += 1
    return no_of_slater_determinants, slater_determinants
    
    
def pairing_matrix_element(no_of_pairs,no_of_levels,slater_det_1,slater_det_2,d,g):
    no_of_same_levels = np.dot(slater_det_1,slater_det_2) #this gives the number levels that are occupied in the 
    #bra-state as well as in the ket-state
    if no_of_same_levels == no_of_pairs: #diagonal matrix elements
        single_particle_energy_factor = sum(np.nonzero(slater_det_1)[0])
        returnvalue = 2*single_particle_energy_factor*d - no_of_pairs*g/2
    elif no_of_same_levels == no_of_pairs - 1: #matrix elements where bra and ket differ in one pair
        returnvalue = -g/2
    else: #Slater-Condon rule
        returnvalue = 0
    return returnvalue
    

def make_pairing_matrix(no_of_pairs,no_of_levels,d,g):
    matrix_size, slater_determinants = make_slater_determinants(no_of_pairs,no_of_levels) #setting up necessary 
    #Slater determinants
    pairing_matrix = np.zeros(shape=(matrix_size,matrix_size)) 
    for i in range(matrix_size):
        for j in range(i,matrix_size):
            pairing_matrix[-i-1,-j-1] = pairing_matrix_element(no_of_pairs,no_of_levels,slater_determinants[i],
                                                               slater_determinants[j],d,g) #matrix elements are 
            #filled up from the bottom such that the lowest states are placed in the beginning
    pairing_matrix = pairing_matrix + pairing_matrix.T - np.diag(pairing_matrix.diagonal()) #makes matrix symmetric 
    #from up-right triangle
    return pairing_matrix


def solve_eigenvalues(no_of_pairs,no_of_levels,d,g_min=-1,g_max=1,no_of_points=21):  
    matrix_size = int(binomial(no_of_levels, no_of_pairs))    
    g_array = np.linspace(g_min,g_max,no_of_points) #values of g for which the problem is solved
    g_mat_eigenvalues = np.zeros(shape=(no_of_points,matrix_size)) #setting up other arrays
    g_mat_eigenvalues_sorted = np.zeros(shape=(no_of_points,matrix_size))
    g_mat_gs_correlation_energy = np.zeros(no_of_points)
    g_mat_eigenvectors = np.zeros(shape=(no_of_points,matrix_size,matrix_size))
    g_mat_ground_state_eigenvector_first_components = np.zeros(no_of_points)

    for x in range(no_of_points): #solving the eigenvalue problems
        g_mat_eigenvalues[x], g_mat_eigenvectors[x] = la.eig(make_pairing_matrix(no_of_pairs,no_of_levels,d,
                                                                                 g_array[x]))
        g_mat_eigenvalues_sorted[x] = np.sort(g_mat_eigenvalues[x])
        g_mat_gs_correlation_energy[x] = g_mat_eigenvalues_sorted[x,0]-(make_pairing_matrix(no_of_pairs,no_of_levels,
                                                                                            d,g_array[x])[0,0]) 

        for i in range(matrix_size): #sort eigenvectors similar as eigenvalues to identify ground state
            if g_mat_eigenvalues_sorted[x,0] == g_mat_eigenvalues[x,i]:
                g_mat_ground_state_eigenvector_first_components[x] = abs(g_mat_eigenvectors[x,0,i])
                break

    #print(g_array,g_mat_gs_correlation_energy)

    plt.figure(figsize=(8,5)) #doing all the plots
    plt.xlabel(r'Interaction strength $g$', fontsize=14)
    plt.ylabel(r'Eigenvalues', fontsize=14)    
    plt.plot(g_array,g_mat_eigenvalues_sorted,"-*")
    #plt.show() 
    plt.savefig("eigenvalues.pdf")

    plt.figure(figsize=(8,5))
    plt.xlabel(r'Interaction strength $g$', fontsize=14)
    plt.ylabel(r'Correlation energy', fontsize=14)    
    plt.plot(g_array,g_mat_gs_correlation_energy,"-*")
    #plt.show() 
    plt.savefig("corr_energy.pdf")

    plt.figure(figsize=(8,5))
    plt.xlabel(r'Interaction strength $g$', fontsize=14)
    plt.ylabel(r'Strength of lowest SD in GS', fontsize=14)
    plt.plot(g_array,g_mat_ground_state_eigenvector_first_components,"-*")
    #plt.show() 
    plt.savefig("gs_strength.pdf")   


def pairing_fci_main(no_of_prt,no_of_states,g):  
    d = 1
    no_of_pairs = int(no_of_prt / 2) #for compatibilty reasons
    no_of_levels = int(no_of_states / 2) #for compatibilty reasons
    
    matrix_size = int(binomial(no_of_levels, no_of_pairs))    
    g_eigenvalues = np.zeros(matrix_size) #setting up other arrays
    g_eigenvalues_sorted = np.zeros(matrix_size)
    g_eigenvectors = np.zeros(shape=(matrix_size,matrix_size))

    #solving the eigenvalue problems
    g_eigenvalues, g_eigenvectors = la.eig(make_pairing_matrix(no_of_pairs,no_of_levels,d,g))
    g_eigenvalues_sorted = np.sort(g_eigenvalues)
    g_gs_correlation_energy = g_eigenvalues_sorted[0]-(make_pairing_matrix(no_of_pairs,no_of_levels,d,g)[0,0]) 

    for i in range(matrix_size): #sort eigenvectors similar as eigenvalues to identify ground state
        if g_eigenvalues_sorted[0] == g_eigenvalues[i]:
            g_ground_state_eigenvector_first_component = abs(g_eigenvectors[0,i])
            break
    
    print("  The FCI result for the correlation energy is ", g_gs_correlation_energy)
    return g_gs_correlation_energy


