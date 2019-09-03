
import numpy as np
from typing import Callable
from abc import ABC, abstractmethod


class PipeWithinPipe:

    @staticmethod
    def get_r_from_xy(x, y, x0=0, y0=0):
        r = np.sqrt(pow(x0 - x, 2) + pow(y0 - y, 2))
        return r

    @abstractmethod
    def get_temperature_r(self, r):
        pass

    def get_temperature_xy(self, x, y, x0=0, y0=0):
        """
        :param x:
        :param y:
        :param x0: origin of radial CSYS - x center of the pipe
        :param y0: origin of radial CSYS - y center of the pipe
        :return:
        """
        r = self.get_r_from_xy(x, y, x0, y0)
        return self.get_temperature_r(r)


class PipeWithinPipeNeumann(PipeWithinPipe):
    def __init__(self, r0, r2, dT_over_dr, T2):
        """
        2D, steady state case:
        Steady heat conduction in a two layered cylinder.
        Temperature at inner (r0) and outer (r2) radius is prescribed.
        r0: float  # inner radius
        r2: float  # outer radius
        J0: float  # temperature gradient (dT/dr) at r = r0. Don't confuse with heat flux, J = k * dT/dr
        T2: float  # temperature for r = r2
        """

        self.r0 = r0
        self.r2 = r2
        self.dT_over_dr = dT_over_dr
        self.T2 = T2

    def get_temperature_r(self, r):
        if r < self.r0:
            return np.nan

        if r > self.r2:
            return np.nan
            # raise ValueError

        T = self.T2 + self.dT_over_dr * self.r0 * (np.log(r / self.r0) - np.log(self.r2 / self.r0))
        return T


class PipeWithinPipeDirichlet(PipeWithinPipe):
    def __init__(self, r0, r1, r2, k1, k2, T0, T2):
        """
        2D, steady state case:
        Steady heat conduction in a two layered cylinder.
        Temperature at inner (r0) and outer (r2) radius is prescribed.
        r0: float  # inner radius
        r1: float  # interface between layers
        r2: float  # outer radius
        k1: float  # inner layer - heat conductivity for r0 < r < r1
        k2: float  # outer layer - heat conductivity for r1 < r < r2
        T0: float  # temperature for r = r0
        T2: float  # temperature for r = r2
        """
        self.r0 = r0
        self.r1 = r1
        self.r2 = r2
        self.k1 = k1
        self.k2 = k2
        self.T0 = T0
        self.T2 = T2

        R1 = np.log(self.r1 / self.r0) / (2 * np.pi * self.k1)  # inner layer - thermal resistance
        R2 = np.log(self.r2 / self.r1) / (2 * np.pi * self.k2)  # outer layer - thermal resistance
        Q = (self.T2 - self.T0) / (R2 + R1)  # heat flux
        self.T1 = Q * R1 + self.T0  # temperature for r = r1 (interface between layers)

    def get_temperature_r(self, r):
        if r < self.r0:
            return self.T0

        if r > self.r2:
            return self.T2
            # raise ValueError

        if r <= self.r1:
            T = (self.T1-self.T0)*np.log(r/self.r0)/np.log(self.r1/self.r0) + self.T0
            return T
        if r <= self.r2:
            T = (self.T2 - self.T1)*np.log(r/self.r1)/np.log(self.r2/self.r1) + self.T1
            return T
