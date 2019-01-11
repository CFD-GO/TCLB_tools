
from SymbolicCollisions.core.cm_symbols import moments_dict
from sympy import exp, pi, integrate, oo
import sympy as sp
import numpy as np
from sympy import diff, ln, sin, pprint
from sympy import Symbol
from sympy.matrices import Matrix
#
# from sympy.interactive.printing import init_printing

from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, \
    F2D, dzeta2D, u2D, \
    rho, w, m00

from SymbolicCollisions.core.printers import round_and_simplify

from joblib import Parallel, delayed
import multiprocessing


class ContinousCMTransforms:
    def __init__(self, dzeta, u, F, rho):
        """
        :param dzeta: direction (x,y,z)
        :param u: velocity (x,y,z)
        :param u: Force (x,y,z)
        :param rho: density (not necessarily m00, for instance in multiphase flows)
        """
        self.dzeta = dzeta
        self.u = u
        self.F = F
        self.rho = rho

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

        cs2 = 1. / 3.
        PI = np.pi
        # PI = sp.pi

        dim = len(self.dzeta)  # number od dimensions
        dzeta_minus_u = self.dzeta - u
        dzeta_u2 = dzeta_minus_u.dot(dzeta_minus_u)

        DF = psi / pow(2 * PI * cs2, dim/2)
        DF *= exp(-dzeta_u2 / (2 * cs2))

        return DF

    def get_hydro_DF(self):
        DF_p = self.get_Maxwellian_DF(psi=(m00 - 1), _u=Matrix([0, 0, 0]))
        DF_gamma = self.get_Maxwellian_DF(psi=1, _u=self.u,)
        return DF_p + DF_gamma

    def get_force_He_hydro_DF(self):
        """
        'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
        """
        cs2 = 1./3.
        eu = self.dzeta.dot(self.F)
        DF_p = self.get_Maxwellian_DF(psi=(m00 - 1), _u=Matrix([0, 0, 0]))

        euF = (self.dzeta - self.u).dot(self.F)
        DF_gamma = self.get_Maxwellian_DF(psi=1, _u=self.u)

        R = -(eu * DF_p + euF * DF_gamma) / (self.rho * cs2)
        R = -R  # `-` sign is skipped to ease code copy-paste ;p
        return R

    def get_force_He_MB(self):
        """
        'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
        Use Maxwellian to calculate equilibria
        """
        cs2 = 1. / 3.
        # cs2 = Symbol('cs2')
        eu_dot_f = (self.dzeta - self.u).dot(self.F)
        result = self.get_Maxwellian_DF() * eu_dot_f / (self.rho * cs2)

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

        cs2 = 1. / 3.
        dim = len(self.dzeta)  # dimension of the space
        w_ = 1. / pow((2 * pi * cs2), dim / 2.)
        w_ *= exp(-e2 / (2 * cs2))
        return w_

    def get_force_Guo(self):
        cs2 = 1. / 3.
        # cs2 = Symbol('cs2')

        eu_terms = self.dzeta - self.u + self.dzeta.dot(self.u)*self.dzeta/cs2
        result = self.get_weight() * self.F.dot(eu_terms) / (self.rho * cs2)
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


def get_mom_vector_from_continuous_def_new(fun, continuous_transformation, moments_order):
    """
    # obviously 2D is faster
    # However 3D works for 2D as well
    :param fun:
    :param continuous_transformation:
    :param moments_order:
    :return:
    """
    # for example: continuous_transformation=get_continuous_cm

    # row = moments_order[0]
    # result = continuous_transformation(row, fun)

    # serial run
    # result = [continuous_transformation(row, fun) for row in moments_order]  # dziala

    # if you experience debugger crashing then run a serial version

    # /pycharm-2018.3.1/helpers/pydev/pydevd.py", line 1487, in dispatch
    #     host = setup['client']
    # TypeError: 'NoneType' object is not subscriptable
    # run in parallel:
    num_cores = multiprocessing.cpu_count()
    result = Parallel(n_jobs=num_cores)(delayed(continuous_transformation)(row, fun) for row in moments_order)

    return Matrix([result])

    #Parallel(n_jobs=2)(delayed(sqrt)(i ** 2) for i in range(10))


def get_mom_vector_from_shift_Mat(fun, Mat):
    pop = Matrix([fun() for i in range(9)])
    # pop = Matrix(9, 1, lambda i,j: i+j)  # column vect
    cm_ = Mat * pop  #for example: Mat=Nraw * Mraw)
    cm_ = round_and_simplify(cm_)
    return Matrix([cm_])

