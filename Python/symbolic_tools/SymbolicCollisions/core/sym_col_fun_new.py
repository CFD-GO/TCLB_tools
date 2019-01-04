
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
from SymbolicCollisions.core.printers import round_and_simplify

from joblib import Parallel, delayed
import multiprocessing

from SymbolicCollisions.core.printers import print_as_vector, print_ccode
import time


def get_continuous_Maxwellian_DF_new(dzeta, u, psi=m00):
    """
    :param dzeta: direction (x,y,z)
    :param u: velocity (x,y,z)
    :param psi: quantity of interest aka scaling function like density
    :return: continuous, local Maxwell-Boltzmann distribution
    'Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments'
    Kannan N. Premnath, Sanjoy Banerjee, 2009
    eq 22
    """

    cs2 = 1. / 3.
    PI = np.pi
    # PI = sp.pi

    D = len(dzeta)  # number od dimensions
    dzeta_u2 = 0
    for dzeta_i, u_i in zip(dzeta, u):
        dzeta_u2 += (dzeta_i-u_i)*(dzeta_i-u_i)

    DF = psi / pow(2 * PI * cs2, D/2)
    DF *= exp(-dzeta_u2 / (2 * cs2))

    return DF


def get_continuous_cm_new(mno, dzeta, u, DF):
    # m, n, o = mno
    # dzeta_x, dzeta_y, dzeta_z = dzeta
    # ux, uy, uz = u

    fun = DF(dzeta, u)
    for dzeta_i, u_i, mno_i in zip(dzeta, u, mno):
        fun *= pow((dzeta_i - u_i), mno_i)
    # fun = DF(dzeta, u) * pow((dzeta_x - ux), m) * pow((dzeta_y - uy), n) * pow((dzeta_z - uz), o)

    lim = [(dim, -oo, oo) for dim in dzeta]
    result = integrate(fun, *lim)

    # result = integrate(fun, (dzeta_x, -oo, oo), (dzeta_y, -oo, oo), (dzeta_z, -oo, oo))
    # result = integrate(fun, (dzeta_x, -oo, oo), (dzeta_y, -oo, oo))
    return round_and_simplify(result)


def get_mom_vector_from_continuous_def_new(fun, continuous_transformation, moments_order):
    # for example: continous_transformation=get_continuous_cm


    # result = get_continuous_cm_new(row, dzeta2D, u2D, get_continuous_Maxwellian_DF_new)  # oczywiscie w 2D liczy szybciej
    # result = get_continuous_cm_new(row, dzeta3D, u3D, get_continuous_Maxwellian_DF_new)
    # result = [get_continuous_cm_new(row, dzeta3D, u3D, get_continuous_Maxwellian_DF_new) for row in moments_order]  # dziala

    # serial run
    # result = [continuous_transformation(row, dzeta3D, u3D, fun) for row in moments_order]  # dziala

    # run in parallel:
    num_cores = multiprocessing.cpu_count()
    result = Parallel(n_jobs=num_cores)(delayed(continuous_transformation)(row, dzeta3D, u3D, fun) for row in moments_order)
    return Matrix([result])

    #Parallel(n_jobs=2)(delayed(sqrt)(i ** 2) for i in range(10))

start = time.process_time()
# time.time() returns a value than can be interpreted as the current date and time.
# time.clock() returns a measure of how much CPU time has been used by the current process.

print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)


cm_eq = get_mom_vector_from_continuous_def_new(get_continuous_Maxwellian_DF_new,
                                               continuous_transformation=get_continuous_cm_new,
                                               moments_order=moments_dict['D2Q9'])

# cm_eq = Matrix([0.9999999999*m00, 0, 0.3333333331*m00, 0.3333333331*m00, 0, 0, 0, 0.1111111112*m00])
print_as_vector(cm_eq, 'cm_eq', regex=True)
print("-----------------------------------------------------------------------------------------------------------")


print(f'\n\n Done in {time.process_time() - start} [s].')
