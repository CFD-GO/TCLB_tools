#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:57:59 2021

@author: grzegorz
"""


import CLB.CLBXMLWriter as CLBXML
import CLB.CLBXMLHandler
import CLB.VTIFile
import os, sys
import numpy as np
from numpy.testing import assert_almost_equal
from numpy.fft import fft2, fftshift, ifft2 # Python DFT
import scipy.optimize as so
import scipy.integrate as sint
import glob
import re
import CLB.VTIFile
import pandas as pd
import math
import matplotlib.pyplot as plt
import time

####
# Show plots in the notebook (don't use it in Python scripts)
# %matplotlib inline 

# WATCH OUT: the implementation of DFT is not 'normalized'
# https://docs.scipy.org/doc/scipy-1.1.0/reference/tutorial/fftpack.html
# https://numpy.org/doc/stable/reference/routines.fft.html#module-numpy.fft

# https://www.unioviedo.es/compnum/labs/PYTHON/lab06_Fourier2D.html


################ SETUP ################



# L = 128
diffusivity0 = 1
tc = 0.01

L = 2
n = 100+1
################ END OF SETUP ################


# Mesh on the square [0,L)x[0,L)

# x = np.linspace(0, L-1, L) + 0.5     # columns (Width)
# y = np.linspace(0, L-1, L) + 0.5     # rows (Height)
# [X,Y] = np.meshgrid(x,y)

# xnorm = (X-0.5) / L * 2 * np.pi
# ynorm = (Y-0.5) / L * 4 * np.pi
# initial_condition = np.exp(np.sin(xnorm)) - 2*np.exp(np.sin(ynorm))


x = np.linspace(0, L, n, endpoint=True) 
y = np.linspace(0, L, n, endpoint=True)
xx, yy = np.meshgrid(x, y)

initial_condition = np.zeros((n, n))


def get_r_from_xy(x, y, x0=0, y0=0):
    r = np.sqrt(pow(x0 - x, 2) + pow(y0 - y, 2))
    return r

for i in range(n):
    for j in range(n):
        r = get_r_from_xy(xx[i][j], yy[i][j], L/2., L/2.)
        if r < 0.25:
            initial_condition[i, j] = 1


plt.figure(figsize=(10, 10))
plt.title(f'initial_condition')
plt.imshow(initial_condition, cmap = 'coolwarm') # coolwarm, gray
fig = plt.gcf()  # get current figure
plt.colorbar()
# plt.close(fig))

F = fft2(initial_condition)  

################ some tests... ################

################ check FFT
                    
Fshift = fftshift(F) # center it in the origin
Fshift = Fshift/(L*L)  # normalize by the area
P = np.abs(Fshift) 
# P = np.real(Fshift)  
# P = np.imag(Fshift)   
 
plt.figure(figsize=(10, 10))                        
plt.imshow(P, extent = [-L,L,-L,L]);
plt.title(f'dft') 
fig = plt.gcf()  # get current figure
plt.colorbar()  
plt.show()
    
    
# ############### check inverse FFT
# yinv = ifft2(F)
# P = np.abs(yinv)    
# P = np.real(yinv)  # TODO: czemu to?
# # P = np.imag(yinv)   

# plt.figure(figsize=(10, 10))     
# plt.imshow(P, cmap = 'gray') # it is wise to check: P-initial_condition
# fig = plt.gcf()  # get current figure
# plt.colorbar()    
# # plt.close(fig)    
    
################ end of tests... ################
   
x2 = np.array([ float(i) if i < n/2. else float(-(n-i)) for i in range(0,n)])
k2, k1 = np.meshgrid(x2, x2)

k1 *= 2.*np.pi/L
k2 *= 2.*np.pi/L
tmp = -diffusivity0 * tc* (k1**2+ k2**2)
decay = np.exp(tmp)

# decay=np.zeros((n,n))
# for m in range(L):
#     # print(f"calculating decay-m: {m}")
#     for n in range(L):
#         # print(f"calculating decay-n: {n}")
#         k1 *= m*2*np.pi/L
#         k2 *= n*2*np.pi/L
#         tmp = -diffusivity0 * tc* (k1**2+ k2**2)
#         decay[m,n] = np.exp(tmp)
    
    
# https://stackoverflow.com/questions/40034993/how-to-get-element-wise-matrix-multiplication-hadamard-product-in-numpy
# decay  = np.ones(decay.shape) # dummy decay
yinv = ifft2(np.multiply(F,decay))

# P = np.abs(yinv)
P = np.real(yinv)   # TODO: czemu to?
# P = np.imag(yinv)   


plt.figure(figsize=(10, 10))    
plt.imshow(P, cmap = 'coolwarm') # it is wise to check: P-initial_condition
plt.title(f'solution after dt') 
fig = plt.gcf()  # get current figure
plt.colorbar()  
plt.show()
# plt.close(fig)

