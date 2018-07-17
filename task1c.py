#this file is basically the same as task1b, only the matrix differs
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

d = 1

garray = np.linspace(-1,1,21)
gmateigenvalues = np.zeros(shape=(21,6))
gmateigenvaluessorted = np.zeros(shape=(21,6))
gmatgscorrelationenergy = np.zeros(21)
gmateigenvectors = np.zeros(shape=(21,6,6))
gmatgroundstateeigenvectorfirstcomponents = np.zeros(21)
for x in range(21):
    g = (x-10)/10
    gmat = -g/2*np.ones((6,6)) #I was too lazy to type the whole matrix and instead I spend way too much time to figure out this "compact" notation
    gmat = gmat-g/2*np.eye(6) 
    gmat = gmat+np.diag([2*d,4*d,6*d,6*d,8*d,10*d])
    for i in range(6):
        gmat[i,5-i] = 0
    #print(gmat)
    gmateigenvalues[x], gmateigenvectors[x] = la.eig(gmat)
    gmateigenvaluessorted[x] = np.sort(gmateigenvalues[x])
    gmatgscorrelationenergy[x] = gmateigenvaluessorted[x,0]-(gmat[0,0]) 
    gmatgroundstateeigenvectorfirstcomponents[x] = abs(gmateigenvectors[x,0,0])  
    #if it is needed to sort eigenvectors accordingly something like this can be used:
    #for i in range(6):
        #if gmateigenvaluessorted[x,0] != gmateigenvalues[x,i]:
            #pass
        #else:
            #print(i)
            #gmatgroundstateeigenvectorfirstcomponents[x] = abs(gmateigenvectors[x,i,0])  
            #break
#print(gmateigenvalues)
#print(gmateigenvaluessorted)

plt.figure()
plt.plot(garray,gmateigenvaluessorted)
plt.title("Eigenvalues as a function of interaction strength")
#plt.show() 
plt.savefig("eigenvalues.png")
plt.figure()
plt.plot(garray,gmatgscorrelationenergy)
plt.title("Correlation energy as a function of interaction strength")
#plt.show() 
plt.savefig("corr_energy.png")
plt.figure()
plt.plot(garray,gmatgroundstateeigenvectorfirstcomponents)
plt.title("Strength of first Slater determinant in ground state")
#plt.show() 
plt.savefig("gs_strength.png")


