
from sympy import exp, pi, integrate, oo
import sympy as sp
import numpy as np

from sympy.matrices import Matrix

from SymbolicCollisions.core.cm_symbols import m00

from SymbolicCollisions.core.printers import round_and_simplify

from joblib import Parallel, delayed
import multiprocessing
from sympy import Symbol
import warnings

class ContinousCMTransforms:
    def __init__(self, dzeta, u, F, rho, cs2=1./3.):
        """
        :param dzeta: direction (x,y,z)
        :param u: velocity (x,y,z)
        :param u: Force (x,y,z)
        :param rho: density (not necessarily m00, for instance in multiphase flows)
        :param cs2: (speed of sound)^2, for isothermal LB cs2=1./3;
                    otherwise  cs2 = Symbol('RT', positive=True)  # positive, negative, real, nonpositive, integer, prime and commutative.

        """
        self.dzeta = dzeta
        self.u = u
        self.F = F
        self.rho = rho
        self.cs2 = cs2

    def get_Maxwellian_DF(self, psi=m00, _u=None):
        """
        :param _u: velocity (x,y,z)
        :param psi: quantity of interest aka scaling function like density
        :return: continuous, local Maxwell-Boltzmann distribution
        'Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments'
        Kannan N. Premnath, Sanjoy Banerjee, 2009
        eq 22
        """
        u = None
        if _u:
            u = _u
        else:
            u = self.u

        PI = sp.pi
        dzeta_minus_u = self.dzeta - u
        dzeta_u2 = dzeta_minus_u.dot(dzeta_minus_u)

        # thank you sympy...
        # hacks:
        # for 2D
        # df = psi / pow(2 * PI * self.cs2, 2/2)
        # LOL: 2/2 gives not simplified result for m22 on d2q9:
        # 1.0 * m00 * (RT * u.y ** 2 - RT ** 1.0 * u.y ** 2 + RT ** 2.0);
        # while 1 does ;p
        # RT ** 2 * m00;

        dim = len(self.dzeta)  # number od dimensions
        # df = psi / pow(2 * PI * self.cs2, dim / 2)  # this is to difficult for sympy :/

        if dim == 2:
            df = psi / (2 * PI * self.cs2)  # 2D version hack
        else:
            df = psi / pow(2 * PI * self.cs2, dim / 2)  # this may be to difficult for sympy :/
            if self.cs2 != 1. / 3.:
                warnings.warn("Sympy may have problem with 3D non isothermal version (cs2=RT) \n "
                              "It also can't simplify it, thus check the raw output", UserWarning)

        df *= exp(-dzeta_u2 / (2 * self.cs2))
        return df

    def get_internal_energy_Maxwellian_DF(self, psi=m00, _u=None):
        df = self.get_Maxwellian_DF(psi=psi, _u=self.u)
        dzeta_minus_u = self.dzeta - self.u
        df = df*0.5*dzeta_minus_u.dot(dzeta_minus_u)
        return df

    def get_total_energy_Maxwellian_DF(self, psi=m00, _u=None):
        df = self.get_Maxwellian_DF(psi=psi, _u=self.u)
        df = df*0.5*self.dzeta.dot(self.dzeta)
        return df

    def get_hydro_DF(self):
        df_p = self.get_Maxwellian_DF(psi=(m00 - 1), _u=Matrix([0, 0, 0]))
        df_gamma = self.get_Maxwellian_DF(psi=1, _u=self.u,)
        return df_p + df_gamma

    def get_force_He_hydro_DF(self):
        """
        'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
        """
        eu = self.dzeta.dot(self.F)
        df_p = self.get_Maxwellian_DF(psi=(m00 - 1), _u=Matrix([0, 0, 0]))

        euf = (self.dzeta - self.u).dot(self.F)
        df_gamma = self.get_Maxwellian_DF(psi=1, _u=self.u)

        R = -(eu * df_p + euf * df_gamma) / (self.rho * self.cs2)
        R = -R  # `-` sign is skipped to ease code copy-paste ;p
        return R

    def get_force_He_MB(self):
        """
        'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
        Use Maxwellian to calculate equilibria
        """
        eu_dot_f = (self.dzeta - self.u).dot(self.F)
        result = self.get_Maxwellian_DF() * eu_dot_f / (self.rho * self.cs2)

        return result

    def get_weight(self):
        """
        PhD Thesis: `The lattice Boltzmann method: Fundamentals and acoustics`
        by Erlend Magnus Viggen
        4.1  The discrete-velocity Boltzmann equation, pp75
        :param i: i-th lattice direction
        :return: returns weight in i-th lattice direction
        """
        e2 = self.dzeta.dot(self.dzeta)

        dim = len(self.dzeta)  # dimension of the space
        w_ = 1. / pow((2 * pi * self.cs2), dim / 2.)
        w_ *= exp(-e2 / (2 * self.cs2))
        return w_

    def get_force_Guo(self):
        eu_terms = self.dzeta - self.u + self.dzeta.dot(self.u)*self.dzeta/self.cs2
        result = self.get_weight() * self.F.dot(eu_terms) / (self.rho * self.cs2)
        return result

    def get_m(self, mno, DF, *args, **kwargs):
        fun = DF(*args, **kwargs)
        for dzeta_i, mno_i in zip(self.dzeta, mno):
            fun *= pow(dzeta_i, mno_i)

        lim = [(dim, -oo, oo) for dim in self.dzeta]
        result = integrate(fun, *lim)
        return round_and_simplify(result)

    def get_cm(self, mno, DF, *args, **kwargs):
        fun = DF(*args, **kwargs)
        for dzeta_i, u_i, mno_i in zip(self.dzeta, self.u, mno):
            fun *= pow((dzeta_i - u_i), mno_i)

        lim = [(dim, -oo, oo) for dim in self.dzeta]
        result = integrate(fun, *lim)

        return round_and_simplify(result)


def get_mom_vector_from_continuous_def(fun, continuous_transformation, moments_order, serial_run=False):
    """
    # obviously 2D is faster
    # However 3D works for 2D as well
    :param fun:
    :param continuous_transformation:
    :param moments_order:
    :param serial_run:
    python debugger may crash in parallel mode.
    Moreover code coverage doesn't work multiprocessing, since the processes are independent beings,
    :return:
    """
    # for example: continuous_transformation=get_continuous_cm

    # row = moments_order[0]
    # result = continuous_transformation(row, fun)

    if serial_run:
        result = [continuous_transformation(row, fun) for row in moments_order]  # serial run
    else:  # run in parallel
        # if you experience debugger crashing then run a serial version
        # /pycharm-2018.3.1/helpers/pydev/pydevd.py", line 1487, in dispatch
        #     host = setup['client']
        # TypeError: 'NoneType' object is not subscriptable

        num_cores = multiprocessing.cpu_count()
        result = Parallel(n_jobs=num_cores)(delayed(continuous_transformation)(row, fun) for row in moments_order)

    return Matrix([result])

    #Parallel(n_jobs=2)(delayed(sqrt)(i ** 2) for i in range(10))
