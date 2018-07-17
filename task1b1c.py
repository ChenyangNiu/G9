import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from matplotlib import interactive
interactive(True)

d = 1

garray = np.linspace(-1,1,21)
gmateigenvalues = np.zeros(shape=(21,2))
#gmatlowereigenvalues = np.zeros(21)
gmateigenvectors = np.zeros(shape=(21,2,2))
gmatgroundstateeigenvectorfirstcomponents = np.zeros(21)
for x in range(21):
    g = (x-10)/10
    gmat = np.array([[-g,-g],[-g,2*d-g]])
    gmateigenvalues[x], gmateigenvectors[x] = la.eig(gmat)
    #gmatlowereigenvalues[x] = gmateigenvalues[x,1]
    gmatgroundstateeigenvectorfirstcomponents[x] = abs(gmateigenvectors[x,1,1])
plt.figure()
plt.plot(garray,gmateigenvalues)
plt.title("Eigenvalues as a function of interaction strength")
#plt.show() 
plt.savefig("eigenvalues.png")
plt.figure()
plt.plot(garray,gmatgroundstateeigenvectorfirstcomponents)
plt.title("Strength of first Slater determinant in ground state")
#plt.show() 
plt.savefig("gs_strength.png")

#print(garray)
#print(gmateigenvalues)
#print(gmateigenvectors)
