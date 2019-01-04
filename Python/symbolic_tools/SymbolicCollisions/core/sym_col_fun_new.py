
from SymbolicCollisions.core.cm_symbols import moments_dict
from sympy import exp, pi, integrate, oo
import sympy as sp
import numpy as np
from sympy import diff, ln, sin, pprint
from sympy import Symbol
from sympy.matrices import Matrix
from sympy.interactive.printing import init_printing
from SymbolicCollisions.core.cm_symbols import ex_D2Q9, ey_D2Q9, w, m00,\
    dzeta3D, u3D, dzeta2D, u2D

from SymbolicCollisions.core.cm_symbols import F3D, F2D, Fx, Fy, F_phi_x, F_phi_y, rho

from SymbolicCollisions.core.printers import round_and_simplify

from joblib import Parallel, delayed
import multiprocessing

from SymbolicCollisions.core.printers import print_as_vector, print_as_vector_old, print_ccode
import time


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
        :param dzeta: direction (x,y,z)
        :param u: velocity (x,y,z)
        :param psi: quantity of interest aka scaling function like density
        :return: continuous, local Maxwell-Boltzmann distribution
        'Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments'
        Kannan N. Premnath, Sanjoy Banerjee, 2009
        eq 22
        """

        if _u:
            u=_u
        else:
            u = self.u


        cs2 = 1. / 3.
        PI = np.pi
        # PI = sp.pi

        dim = len(self.dzeta)  # number od dimensions
        dzeta_u2 = 0
        for dzeta_i, u_i in zip(self.dzeta, u):
            dzeta_u2 += (dzeta_i-u_i)*(dzeta_i-u_i)

        DF = psi / pow(2 * PI * cs2, dim/2)
        DF *= exp(-dzeta_u2 / (2 * cs2))

        return DF

    def get_hydro_DF(self):
        # DF_p = get_continuous_Maxwellian_DF(dzeta=(dzeta_x, dzeta_y), psi=(m00 - 1), u=(0, 0))
        # DF_gamma = get_continuous_Maxwellian_DF(dzeta=(dzeta_x, dzeta_y), psi=1, u=(ux, uy))
        DF_p = self.get_Maxwellian_DF(psi=(m00 - 1), _u=(0, 0, 0))
        DF_gamma = self.get_Maxwellian_DF(psi=1, _u=self.u,)
        return DF_p + DF_gamma

    def get_force_He_hydro_DF_2D(self):
        """
        'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
        """
        # cs2 = 1. / 3.
        # eu = dzeta[0] * Fx + dzeta[1] * Fy
        # dzeta_2 = dzeta[0] * dzeta[0] + dzeta[1] * dzeta[1]
        # DF_p = (m00 - 1) / (2 * pi * cs2) * exp(-dzeta_2 / (2 * cs2))
        # # DF_p = get_continuous_Maxwellian_DF(dzeta=dzeta, psi=(m00-1), u=(0, 0))
        #
        # euF = (dzeta[0] - ux) * Fx + (dzeta[1] - uy) * Fy
        # dzeta_u_2 = (dzeta[0] - ux) * (dzeta[0] - ux) + (dzeta[1] - uy) * (dzeta[1] - uy)
        # DF_gamma = 1 / (2 * pi * cs2) * exp(-dzeta_u_2 / (2 * cs2))
        # # DF_gamma = get_continuous_Maxwellian_DF(dzeta=dzeta, psi=1, u=(ux, uy))
        #
        # R = -(eu * DF_p + euF * DF_gamma) / (rho * cs2)
        # R = -R  # `-` sign is skipped to ease code copy-paste ;p
        # return R

        cs2 = 1. / 3.
        eu = self.dzeta[0] * self.F[0] + self.dzeta[1] * self.F[1]
        dzeta_2 = self.dzeta[0] * self.dzeta[0] + self.dzeta[1] * self.dzeta[1]
        DF_p = (m00 - 1) / (2 * pi * cs2) * exp(-dzeta_2 / (2 * cs2))
        # DF_p = get_continuous_Maxwellian_DF(dzeta=dzeta, psi=(m00-1), u=(0, 0))

        euF = (self.dzeta[0] - self.u[0]) * self.F[0] + (self.dzeta[1] - self.u[1]) * self.F[1]
        dzeta_u_2 = (self.dzeta[0] - self.u[0]) * (self.dzeta[0] - self.u[1]) + (self.dzeta[1] - self.u[1]) * (self.dzeta[1] - self.u[1])
        DF_gamma = 1 / (2 * pi * cs2) * exp(-dzeta_u_2 / (2 * cs2))
        # DF_gamma = get_continuous_Maxwellian_DF(dzeta=dzeta, psi=1, u=(ux, uy))

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

        euF = 0
        for dzeta_i, u_i, F_i in zip(self.dzeta, self.u, self.F):
            euF += (dzeta_i - u_i) * F_i

        R = self.get_Maxwellian_DF() * euF / (self.rho * cs2)

        return R

    def get_weight(self):
        """
        PhD Thesis: `The lattice Boltzmann method: Fundamentals and acoustics`
        by Erlend Magnus Viggen
        4.1  The discrete-velocity Boltzmann equation, pp75
        :param i: i-th lattice direction
        :return: returns weight in i-th lattice direction
        """

        e2 = 0
        for dzeta_i in self.dzeta:
            e2 += dzeta_i * dzeta_i

        cs2 = 1. / 3.
        dim = len(self.dzeta)  # dimension of the space
        w_ = 1. / pow((2 * pi * cs2), dim / 2.)
        w_ *= exp(-e2 / (2 * cs2))
        return w_


    def get_force_Guo_2D(self):
        cs2 = 1. / 3.
        # cs2 = Symbol('cs2')
        # extended version with second order terms

        # TODO write 3D version
        # D = len(self.dzeta)  # number od dimensions
        # temp = [None] * D
        # for i in range(D):
        #     temp[i] = ...

        temp_x = self.dzeta[0] - self.u[0] + (self.dzeta[0] * self.u[0] + self.dzeta[1] * self.u[1]) * self.dzeta[0] / cs2
        temp_y = self.dzeta[1] - self.u[1] + (self.dzeta[0] * self.u[0] + self.dzeta[1] * self.u[1]) * self.dzeta[1] / cs2

        result = self.get_weight() * (temp_x * self.F[0] + temp_y * self.F[1]) / (self.rho * cs2)

        return result

    def get_m(self, mno, DF):
        fun = DF()
        for dzeta_i, mno_i in zip(self.dzeta, mno):
            fun *= pow(dzeta_i, mno_i)

        lim = [(dim, -oo, oo) for dim in self.dzeta]
        result = integrate(fun, *lim)
        return round_and_simplify(result)

    def get_cm(self, mno, DF):
        fun = DF()
        for dzeta_i, u_i, mno_i in zip(self.dzeta, self.u, mno):
            fun *= pow((dzeta_i - u_i), mno_i)

        lim = [(dim, -oo, oo) for dim in self.dzeta]
        result = integrate(fun, *lim)

        return round_and_simplify(result)


def get_mom_vector_from_continuous_def_new(fun, continuous_transformation, moments_order):
    # for example: continous_transformation=get_continuous_cm

    # row = moments_order[0]
    # result = continuous_transformation(row, fun)  # oczywiscie w 2D liczy szybciej
    # result = continuous_transformation(row, fun)
    # result = [continuous_transformation(row, fun) for row in moments_order]  # ale 3d dziala tez dla 2D
    #
    # test=continuous_transformation(row, fun)

    # serial run
    # result = [continuous_transformation(row, fun) for row in moments_order]  # dziala

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


lattice = 'D2Q9'
start = time.process_time()
# time.time() returns a value than can be interpreted as the current date and time.
# time.clock() returns a measure of how much CPU time has been used by the current process.


cm_i = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
# cm_i = ContinousCMTransforms(dzeta2D, u2D, F2D, rho)
# row = moments_dict['D2Q9'][0]
# test= cm_i.get_cm(row, cm_i.get_Maxwellian_DF)
# test=Matrix([test])
# print_as_vector(test, 'test', regex=True)

print("-----------------------------------------------------------------------------------------------------------")
# print('\n//population_eq -> cm_eq - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def_new(cm_i.get_Maxwellian_DF,
#                                                continuous_transformation=cm_i.get_cm,
#                                                moments_order=moments_dict[lattice])


# TODO: get_cm wychodzi jak get_m, popraw przekazywanie argumentow...
cm_eq = get_mom_vector_from_continuous_def_new(cm_i.get_hydro_DF,
                                               continuous_transformation=cm_i.get_cm,
                                               moments_order=moments_dict[lattice])

#
print_as_vector(cm_eq, 'cm_eq', regex=True)
# print_as_vector(cm_eq, 'cm_eq_no_regex', regex=False)

# print('\n//Force -> Force_cm - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = forceM(rho,u,x,y) *(x-ux)^m (y-uy)^n ')
# # F_cm = get_mom_vector_from_continuous_def_new(cm_i.get_force_He_MB,
# F_cm = get_mom_vector_from_continuous_def_new(cm_i.get_force_He_hydro_DF_2D,
#                                               continuous_transformation=cm_i.get_cm,
#                                               moments_order=moments_dict[lattice])
# print_as_vector(F_cm, 'F_cm', regex=True)
#
# print(f'\n\n Done in {time.process_time() - start} [s].')
#






# print('\n//Force -> Force_cm - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = forceM(rho,u,x,y) *(x-ux)^m (y-uy)^n ')
# F_cm_Guo_extended = get_mom_vector_from_discrete_def(get_discrete_force_Guo, discrete_transform=get_discrete_cm)
# print_as_vector(F_cm_Guo_extended, 'F_cm', regex=True)


# print('\n//N*M*F_He_continous ')
# from SymbolicCollisions.core.cm_symbols import \
#     ex_D3Q19 as ex, \
#     ey_D3Q19 as ey, \
#     ez_D3Q19 as ez

# from SymbolicCollisions.core.cm_symbols import \
#     ex_D2Q9 as ex, \
#     ey_D2Q9 as ey
# ez = None
# from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix
#
# Mraw = get_raw_moments_matrix(ex, ey, ez)
# Nraw = get_shift_matrix(Mraw.inv(), ex, ey, ez)
#
# TO TYLKO DLA DYSKRETNEJ WERSJI!!!
# NMF_cm_He_original = get_mom_vector_from_shift_Mat(get_continuous_force_He_MB, Mat=Nraw * Mraw)
# print_as_vector(NMF_cm_He_original, 'F_cm', regex=True)  # produces looong expressions
#