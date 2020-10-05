"""
The analytical solution is from
'Boundary condition at a two-phase interface in the lattice Boltzmann method for the convection-diffusion equation'
Hiroaki Yoshida,Takayuki Kobayashi,Hidemitsu Hayashi,Tomoyuki Kinjo, Hitoshi Washizu and Kenji Fukuzawa, 2014


However, there are nicer plots in
'Phase interface effects in the total enthalpy-based lattice Boltzmann model for solid–liquid phase change'
Rongzong Huang, Huiying Wu, 2015
"""

import numpy as np
from scipy import special

class Solver:
    def __init__(self, T_left, T_right, k_left, k_right, cp_left, cp_right):
        """
        Follows Fig.1. from 'Phase interface effects in the total enthalpy-based lattice Boltzmann model for solid–liquid phase change'
        The transient temperature profile across two layer composite slab (oriented vertically) are to be calculated.
        The density is assumed to be uniform.
        :param T_left: Temperature at the surface of the left layer [K]
        :param T_right: Temperature at the surface of the right layer [K]
        :param k_left: thermal conductivity of the left layer [W/m2K]
        :param k_right: thermal conductivity of the right layer [W/m2K]
        :param cp_left: specyfic heat capacity of the left layer
        :param cp_right: specyfic heat capacity of the right layer [J/(K kg)]
        """

        self.T_left = T_left
        self.T_right = T_right
        self.k_left = k_left
        self.k_right = k_right
        self.cp_left = cp_left
        self.cp_right = cp_right

    def calc_transient_T_profile(self, t, x):
        """
        'Boundary condition at a two-phase interface in the lattice Boltzmann method for the convection-diffusion equation'
        Hiroaki Yoshida,Takayuki Kobayashi,Hidemitsu Hayashi,Tomoyuki Kinjo, Hitoshi Washizu and Kenji Fukuzawa, 2014
        The time step is (KA/λAL2)?t = 0.125 ×
        :param t: time
        :param x: 0 is located between stabs
        :return: Temperature(t,x)
        """

        # back to R Huang notation
        # coeff = np.sqrt((self.k_left * self.cp_left) / (self.k_right * self.cp_right))
        coeff = np.sqrt((self.k_left * self.cp_left) / (self.k_right * self.cp_right))
        if x < 0:
            # B - po lewej - eq 45
            err_arg = -x/(2 * np.sqrt(t * self.k_left / self.cp_left))
            err_result = special.erfc(err_arg)  # special.erfc(err_arg) =  1 - special.erf(err_arg)
        else:
            # A - po prawej - eq 44
            err_arg = x/(2 * np.sqrt(t * self.k_right / self.cp_right))
            err_result = 1 + coeff * special.erf(err_arg)

        T = self.T_left + err_result * (self.T_right - self.T_left) / (1 + coeff)
        return T

    def calc_transient_T_profile_RHuang(self, t, x):
        """
        'Phase interface effects in the total enthalpy-based lattice Boltzmann model for solid–liquid phase change'
        Rongzong Huang, Huiying Wu, 2015
        :param t: time
        :param x: 0 is located between stabs
        :return: Temperature(t,x)
        """

        Warning("This is buggy for R_cp")
        if x < 0:
            err_arg = -x
            err_arg /= (2 * np.sqrt(t*self.k_left / self.cp_left))
            err_result = special.erfc(err_arg)

            coeff = (self.T_left - self.T_right)*np.sqrt(self.k_right*self.cp_right)
            coeff /= (np.sqrt(self.k_left * self.cp_right) + np.sqrt(self.k_right * self.cp_right))

            T = self.T_left - coeff*err_result
            return T

        else:
            err_arg = x
            err_arg /= (2 * np.sqrt(t*self.k_right / self.cp_right))
            err_result = special.erfc(err_arg)

            coeff = (self.T_left - self.T_right) * np.sqrt(self.k_left*self.cp_left)
            coeff /= (np.sqrt(self.k_left * self.cp_right) + np.sqrt(self.k_right * self.cp_right))

            T = self.T_right + coeff * err_result
            return T