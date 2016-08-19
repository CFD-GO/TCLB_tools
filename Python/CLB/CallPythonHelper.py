# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:23:48 2016

@author: mdzik
"""

import numpy as np

class CallPythonHelper:
    def __init__(self, *args):
        self.offsets, self.time = args[:2]
        self.data = args[2:]

    def getVector(self,idx):
        V = np.asarray(self.data[idx])
        shape = V.shape
        return np.asarray(self.data[idx]).reshape((shape[1],shape[0],3))

    def getScalar(self,idx):
        V = np.asarray(self.data[idx])
        shape = V.shape
        return np.asarray(self.data[idx]).reshape((shape[1],shape[0]))

    def getXY(self,scal_idx=0):
        V = np.asarray(self.data[scal_idx])
        shape = V.shape

        X = np.zeros([shape[1],shape[0]])
        Y = np.zeros([shape[1],shape[0]])

        for x in range(shape[0]):
            X[:,x] = x

        for y in range(shape[1]):
            Y[y,:] = y
        return X,Y
    def dump(self, fname):
        ddump = list()
        for i,d in enumerate(self.data):
            shape =  np.asarray(d).shape
            if shape[2] == 1:
                ddump.append( self.getScalar(i) )
            else:
                ddump.append( self.getVector(i) )

        np.savez(fname, ddump)   