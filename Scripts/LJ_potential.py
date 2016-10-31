#!/usr/bin/python

"""
A script to play around with the lennard jones potential
e (epsilon) = well depth
s (sigma) = radius at which the potential is zero
"""
import numpy as np
import matplotlib.pyplot as plt


def vlj(e1, s1, e2, s2, r):
    # define equation for calculation of lennard jones potential
    s = (s1*s2)**.5
    e = (e1*e2)**.5
    V = []
    for i in range(0, len(r)):
        V.append(4*e*((s/r[i])**12 - (s/r[i])**6))
    return V, e, s

e_ref = .978638
s_ref = .3401

r = np.linspace(s_ref - 0.02, 1, 100)

e2 = 2
s2 = .3401


V, e, s = vlj(e2, s2, e_ref, s_ref, r)

rm = s*(2**(1.0/6.0))  # radius at which minimum potential well depth is located

plt.plot(r, V)
plt.plot(rm, -e, 'o')
plt.xlabel('Distance Apart (nm)')
plt.ylabel('Potential (kJ/mol)')
plt.show()