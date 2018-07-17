#!/usr/bin/python
from sympy import *
from pylab import *
import matplotlib.pyplot as plt


ga = linspace(-1,1,20)
e1 = []
e2 = []


for g_val in ga:
    H1 = matrix([[-g_val , -g_val],
                 [-g_val    , 2-g_val]])

    u1, v1 = linalg.eig(H1)
    # e1.append(min(u1))
    c1 = sorted(u1)
    e1.append(c1[0])
    e2.append(c1[1])

    
exact = e1 - (-ga)

plt.axis([-1,1,-2,4])
plt.xlabel(r'Interaction strength, $g$', fontsize=16)
plt.ylabel(r'eigenvalues', fontsize=16)

plt.plot(ga, e1,'b-*',linewidth = 2.0, label = 'e1', color = "red")
plt.plot(ga, -ga,'b-*',linewidth = 2.0, label = '-g',  color = "green")
plt.plot(ga, exact,'b-*',linewidth = 2.0, label = 'corr', color = "blue")

plt.plot(ga, e2,'b-*',linewidth = 2.0, label = 'e2', color = "black", linestyle = "dashed")


plt.legend()
plt.savefig('part1b_niu.pdf', format='pdf')
plt.show()

