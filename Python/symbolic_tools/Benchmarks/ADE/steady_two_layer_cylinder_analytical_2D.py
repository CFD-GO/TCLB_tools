
import numpy as np


class PipeWithinPipe:
    def __init__(self, r0, r1, r2, k1, k2, T0, T2):
        """
        2D, steady state case:
        Steady heat conduction in a two layered cylinder.
        Temperature at inner (r0) and outer (r2) radius is prescribed.

        :param r0: inner radius
        :param r1: interface between layers
        :param r2: outer radius
        :param k1: inner layer - heat conductivity for r0 < r < r1
        :param k2: outer layer - heat conductivity for r1 < r < r2
        :param T0: temperature for r = r0
        :param T2: temperature for r = r2
        """

        self.r0 = r0
        self.r1 = r1
        self.r2 = r2
        self.k1 = k1
        self.k2 = k2
        self.T0 = T0
        self.T2 = T2

        R1 = np.log(r1/r0)/(2*np.pi*k1)  # inner layer - thermal resistance
        R2 = np.log(r2/r1)/(2*np.pi*k2)  # outer layer - thermal resistance
        Q = (T2-T0)/(R2+R1)  # heat flux
        self.T1 = Q*R1 + T0  # temperature for r = r1 (interface between layers)
        print(self.T1)

    def get_temperature(self, r):

        if r < self.r0 or r > self.r2:
            raise ValueError

        if r <= self.r1:
            T = (self.T1-self.T0)*np.log(r/self.r0)/np.log(self.r1/self.r0) + self.T0
            return T
        if r <= self.r2:
            T = (self.T2 - self.T1)*np.log(r/self.r1)/np.log(self.r2/self.r1) + self.T1
            return T
