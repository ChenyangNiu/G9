#this file computes the pairing matrix for a given number of pairs and levels and solves the eigenvalue problem
import numpy as np
from numpy import linalg as la
from sympy.utilities.iterables import multiset_permutations as mup
import matplotlib.pyplot as plt
import math
#import itertools as it


def binomial(a,b):
    return math.factorial(a)/(math.factorial(a-b)*math.factorial(b))


def makeslaterdeterminants(noofpairs,nooflevels):
    noofslaterdeterminants = int(binomial(nooflevels, noofpairs))
    creationarray = np.concatenate((np.ones(noofpairs), np.zeros(nooflevels-noofpairs)))
    #print(noofslaterdeterminants)
    slaterdeterminants = np.zeros(shape=(noofslaterdeterminants,nooflevels))
    i = 0
    for p in mup(creationarray):
        slaterdeterminants[i] = p
        i += 1
    #print(slaterdeterminants)
    return noofslaterdeterminants, slaterdeterminants
    
    
def pairingmatrixelement(noofpairs,nooflevels,slaterdet1,slaterdet2,d,g):
    noofsamelevels = np.dot(slaterdet1,slaterdet2)
    if noofsamelevels == noofpairs:
        singleparticleenergyfactor = sum(np.nonzero(slaterdet1)[0])
        returnvalue = 2*singleparticleenergyfactor*d - noofpairs*g/2
    elif noofsamelevels == noofpairs - 1:
        returnvalue = -g/2
    else:
        returnvalue = 0
    #print(returnvalue)
    return returnvalue
    

def makepairingmatrix(noofpairs,nooflevels,d,g):
    matrixsize, slaterdeterminants = makeslaterdeterminants(noofpairs,nooflevels)
    pairingmatrix = np.zeros(shape=(matrixsize,matrixsize)) 
    for i in range(matrixsize):
        for j in range(i,matrixsize):
            pairingmatrix[-i-1,-j-1] = pairingmatrixelement(noofpairs,nooflevels,slaterdeterminants[i],slaterdeterminants[j],d,g)
    #print(pairingmatrix)
    pairingmatrix = pairingmatrix + pairingmatrix.T - np.diag(pairingmatrix.diagonal()) #makes matrix symmetric from up-right triangle
    #print(pairingmatrix)
    return pairingmatrix


def solveeigenvalues(noofpairs,nooflevels,d):
    matrixsize = int(binomial(nooflevels, noofpairs))    
    garray = np.linspace(-1,1,21)
    gmateigenvalues = np.zeros(shape=(21,matrixsize))
    gmateigenvaluessorted = np.zeros(shape=(21,matrixsize))
    gmatgscorrelationenergy = np.zeros(21)
    gmateigenvectors = np.zeros(shape=(21,matrixsize,matrixsize))
    gmatgroundstateeigenvectorfirstcomponents = np.zeros(21)
    for x in range(21):
        g = (x-10)/10
        gmateigenvalues[x], gmateigenvectors[x] = la.eig(makepairingmatrix(noofpairs,nooflevels,d,g))
        gmateigenvaluessorted[x] = np.sort(gmateigenvalues[x])
        gmatgscorrelationenergy[x] = gmateigenvaluessorted[x,0]-(makepairingmatrix(noofpairs,nooflevels,d,g)[0,0]) 
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

solveeigenvalues(2,4,1)
 


