#!/usr/bin/python
from sympy import *
from pylab import *
import matplotlib.pyplot as plt


ga = linspace(-1,1,20)
e1 = []
e2 = []
e3 = []
e4 = []
e5 = []
e6 = []

for g_val in ga:
    H1 = matrix([[2-g_val , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,     0],
                 [-g_val/2.,   4-g_val,  -g_val/2., -g_val/2.,    0., -g_val/2.],
                 [-g_val/2., -g_val/2.,    6-g_val,     0, -g_val/2., -g_val/2.],
                 [-g_val/2., -g_val/2.,      0,   6-g_val, -g_val/2., -g_val/2.],
                 [-g_val/2.,     0,  -g_val/2., -g_val/2.,   8-g_val, -g_val/2.],
                 [0    , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,  10-g_val]])

    u1, v1 = linalg.eig(H1)
    # e1.append(min(u1))
    c1 = sorted(u1)
    e1.append(c1[0])
    e2.append(c1[1])
    e3.append(c1[2])
    e4.append(c1[3])
    e5.append(c1[4])
    e6.append(c1[5])
    
exact = e1 - (2-ga)

plt.axis([-1,1,-2,12])
plt.xlabel(r'Interaction strength, $g$', fontsize=16)
plt.ylabel(r'eigenvalues', fontsize=16)

plt.plot(ga, e1,'b-*',linewidth = 2.0, label = 'e1', color = "red")
plt.plot(ga, 2-ga,'b-*',linewidth = 2.0, label = '2-g',  color = "green")
plt.plot(ga, exact,'b-*',linewidth = 2.0, label = 'corr', color = "blue")

plt.plot(ga, e2,'b-*',linewidth = 2.0, label = 'e2', color = "red", linestyle = "dashed")
plt.plot(ga, e3,'b-*',linewidth = 2.0, label = 'e3', color = "blue", linestyle = "dashed")
plt.plot(ga, e4,'b-*',linewidth = 1.0, label = 'e4', color = "green")
plt.plot(ga, e5,'b-*',linewidth = 2.0, label = 'e5', color = "black", linestyle = "dashed")
plt.plot(ga, e6,'b-*',linewidth = 2.0, label = 'e6', color = "yellow", linestyle = "dashed")

plt.legend()
plt.savefig('part1c_niu.pdf', format='pdf')
plt.show()
