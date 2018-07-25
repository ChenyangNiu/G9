import numpy as np
import matplotlib.pyplot as plt
import pairing_FCI as fci #task 1
import pairing_CCD as ccd #task 2
import sys
sys.setrecursionlimit(2000) #this allows to change the python-build-in recursion limit

#ccd.pairing_ccd_main(4, 12, -.1, .00001, 100)
#fci.pairing_fci_main(4, 12, .2)

def compare_methods(no_of_prt, no_of_states, g_min, g_max, no_of_g_points, ccd_accuracy, ccd_maxit):    
    g_array = np.linspace(g_min, g_max, no_of_g_points)  
    fci_E_C = np.zeros(no_of_g_points)  
    ccd_E_C = np.zeros(no_of_g_points)
    
    for x in range(no_of_g_points):
        fci_E_C[x] = fci.pairing_fci_main(no_of_prt, no_of_states, g_array[x])    
        ccd_E_C[x] = ccd.pairing_ccd_main(no_of_prt, no_of_states, g_array[x], ccd_accuracy, ccd_maxit)

    plt.figure() #figsize=(8,5)
    plt.xlabel(r'Interaction strength $g$', fontsize=14)
    plt.ylabel(r'Correlation energy', fontsize=14)    
    plt.plot(g_array,fci_E_C,"r-*", label="FCI")  
    plt.plot(g_array,ccd_E_C,"b-*", label="CCD")
    plt.legend(prop={'size': 12})
    #plt.show() 
    plt.savefig("fci_ccd_corr_energy.pdf")


compare_methods(4, 8, -1, 2.2, 33, .00001, 1000)