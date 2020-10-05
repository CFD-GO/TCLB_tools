#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 10:39:54 2018

@author: mdzikowski
"""


class CLBInput:
    def __init__(self):
        pass
    
    


import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
from scipy import ndimage

data = misc.imread('/home/mdzikowski/projekty/ogniwa/ak9_ak10_ak10_800_nh_recon_bh0.8_export.to-byte035.bmp')


data2 = np.where(data == 255, 0 ,1)[200:400,100:300]


NN = data2.size

filter_size = 50

labeled_array, num_features = ndimage.label(data2)


plt.matshow(labeled_array)
plt.colorbar()
plt.matshow(data2)


plt.figure()

sizes = [ NN - np.count_nonzero(labeled_array - l) for l in range(num_features) ]
plt.plot(np.sort(sizes), 'o')
plt.grid(which='both')


data3 = labeled_array.copy()
for l,s in zip( range(num_features), sizes) :
    if s < filter_size:
        data3 = np.where(data3 == l, 0, data3)

np.savez('/tmp/cache.npz', data3=data3)

data3 = np.load('/tmp/cache.npz')['data3']

plt.matshow(data3)

data4 = np.where(data3 > 0, 1 ,0)
plt.matshow(data4)

#init
data5 = np.zeros_like(data4)

data5[:,0] = np.where(data4[:,0] == 0, 1, data5[:,0])


s = data5[0,0] 
cl = 1
mark = list()
mark.append(0)
lmax = 0
markmax = 0
l = 0
for i in range(1,data5[:,0].size):
    if data5[i,0] != data5[i-1,0]:
        if l > lmax:
            lmax = l
            markmax = cl
        l = 0
        cl = cl + 1
    if data5[i,0] == 0:
        l = l + 1
        mark.append(cl)
    else:
        mark.append(0)

plt.figure()
print markmax
print lmax

data5[:,0] = np.where(np.array(mark) == markmax, 1, 0)

print data5.shape

plt.matshow(data5)
plt.colorbar()
plt.show()
ex,ey = np.meshgrid([-1,0,1],[-1,0,1])
idx = np.where(ex**2+ey**2 == 1)
ex = ex[idx]
ey = ey[idx]
cov = 0
while cov - np.sum(data5) != 0:
#for i in range(2):
    print(cov)
    cov = np.sum(data5)
    tmp = np.copy(data5)
    for x_, y_ in zip( ex.ravel(), ey.ravel() ):
        tmp = tmp + np.roll(np.roll(data5,axis=0, shift=x_),axis=1, shift=y_)
        
    data5 = np.where(data4 == 0, tmp,0)    
    data5 = np.where(data5 > 0, 1,0)
    
    
fout = file('/tmp/voxelized.txt', 'w')
foutT = file('/tmp/voxelized-t.txt', 'w')

np.savetxt('/tmp/voxelized.txt', data5, delimiter='\n', fmt='%d')
np.savetxt('/tmp/voxelized-t.txt', data5.T, delimiter='\n', fmt='%d')

plt.matshow(data5)




data = np.loadtxt('/tmp/voxelized-t.txt')
print data.shape
Nx = 200
Ny = data.shape[0] / Nx
print Nx*Ny
data = data.reshape([Nx,Ny])
#plt.matshow(data)
#plt.show()

subdata = data
data2 = ndimage.zoom(subdata, 2, order=0)
plt.matshow( data2 )

ex,ey = np.meshgrid([-1,0,1],[-1,0,1])

data3 = np.copy(data2)

diff = np.zeros_like(data2)

for i in range(4):
    tmp = np.copy(data3)
    for x_, y_ in zip( ex.ravel(), ey.ravel() ):
       tmp = tmp + np.roll(np.roll(data3,axis=0, shift=x_),axis=1, shift=y_)
    #diff = (data3 - tmp)
    data3 = np.where(tmp > 3, 0, 1)
    
    plt.matshow(data3)
    
print data3.shape
np.savetxt('/tmp/voxelized-resampled.txt', np.array(np.where(data3==0,int(1),int(0)).T,dtype=int), delimiter='\n', fmt='%d')
plt.show()
    