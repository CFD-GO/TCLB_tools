#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 13:06:24 2020

@author: grzegorz
"""

import os
import numpy as np
import scipy.optimize as so
import scipy.integrate as sint
import glob
import re
import CLB.VTIFile
import pandas as pd
import matplotlib.pyplot as plt


"""
To fix ideas:
tau = 3v + 1/2
omega = 1/tau

0 < v < 1/6
1/2 < tau < 1
2 > omega > 1

see chapter 10.7.2, eq 10.48, p429 from 'The Lattice Boltzmann Method: Principles and Practice'
by T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
There are certain values of magic_parameter that show distinctive properties:
• magic_parameter 1./12 = 0.08(3) cancels the third-order spatial error, leading to optimal results for pure advection problems.
• magic_parameter 1./6 = 0.1(6) cancels the fourth-order spatial error, providing the most accurate results for the pure diffusion equation.
• magic_parameter 3./16 = 0.1875 results in the boundary wall location implemented via bounce-back for the Poiseuille flow exactly in the middle between horizontal walls and fluid nodes.
• magic_parameter 1./4 = 0.25 provides the most stable simulations.

"""
omega_odd = np.linspace(1,2,100) 
magic_parameter = 1/4.


def get_another_omega(omega_in, magic_parameter):
    omega_out  = 2*(2-omega_in)/(omega_in*(4-magic_parameter-1) + 2)
    return omega_out
    
def get_tau_from_omega(omega):
    return 1./omega

fig = plt.figure()


plt.plot(omega_odd,get_another_omega(omega_odd, magic_parameter=1./4),  label=f'$ \lambda = 1/4$')
plt.plot(omega_odd,get_another_omega(omega_odd, magic_parameter=3./16), label=f'$ \lambda = 3/16$')
plt.plot(omega_odd,get_another_omega(omega_odd, magic_parameter=1./6),  label=f'$ \lambda = 1/6$')
plt.plot(omega_odd,get_another_omega(omega_odd, magic_parameter=1./12), label=f'$ \lambda = 1/12$')


plt.grid(which='both')
plt.xlabel('$\omega_{odd} $')
plt.ylabel('$\omega_{even} $')

plt.legend()

plt.savefig(f'magic_plot.png', dpi=200)

plt.show()
# plt.close(fig)
print("Done.")