"""
M - distributions to raw moment transformation matrix
N - raw moments to central moments transformation matrix

based on:
'Modeling incompressible thermal flows using a central-moment-based lattice Boltzmann method'
Linlin Fei, Kai Hong Luo, Chuandong Lin, Qing Li
2017
"""


from sympy import exp, pi, integrate, oo
from sympy import Symbol
from sympy.matrices import Matrix
from sympy.interactive.printing import init_printing
from SymbolicCollisions.core.cm_symbols import ux, uy, w, m00, \
    Fx, Fy, F_phi_x, F_phi_y, rho, dzeta_x, dzeta_y, \
    NrawD2Q9, Mraw_D2Q9, M_ortho_GS, \
    F2D, F3D, dzeta2D, dzeta3D, u2D, u3D, \
    ex_D2Q9 as ex, \
    ey_D2Q9 as ey


from SymbolicCollisions.core.printers import round_and_simplify

init_printing(use_unicode=False, wrap_line=False, no_global=True)


def get_DF(q=9, print_symbol='default_symbol2'):
    symbols_ = [Symbol("%s[%d]" % (print_symbol, i)) for i in range(0, q)]
    return Matrix(symbols_)


def get_m00(q, print_symbol='default_symbol3'):
    m00_ = Symbol("%s[%d]" % (print_symbol, 0))

    for i in range(1, q):
        m00_ += Symbol("%s[%d]" % (print_symbol, i))

    return m00_


def get_e():
    symbols_ = [Matrix([ex[i], ey[i]]) for i in range(9)]
    return Matrix([symbols_])


def get_gamma(i):
    cs2 = 1. / 3.
    # cs2 = Symbol('cs2')
    eu = ex[i] * ux + ey[i] * uy
    u2 = ux * ux + uy * uy
    gamma = w[i] * (1 + eu / cs2 + eu * eu / (2 * cs2 * cs2) - u2 / (2 * cs2))
    return gamma


def get_discrete_EDF_hydro(i):
    gamma = get_gamma(i)
    g = m00 * w[i] + gamma - w[i]
    return g


def get_discrete_force_He(i):
    """
    'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
    """
    cs2 = 1. / 3.
    # cs2 = Symbol('cs2')
    euF = (ex[i] - ux) * Fx + (ey[i] - uy) * Fy
    pop_eq = m00 * get_gamma(i)
    R = pop_eq * euF / (rho * cs2)
    return R


def get_discrete_force_He_hydro_eq_experimental(i):
    """
    'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
    version for 'Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios' A. Fakhari et. al., 2017
    """
    cs2 = 1. / 3.
    # cs2 = Symbol('cs2')
    euF = (ex[i] - ux) * Fx + (ey[i] - uy) * Fy
    pop_eq = get_discrete_EDF_hydro(i)
    # R = pop_eq*euF/(p_star*cs2)
    R = pop_eq * euF / (rho * cs2)
    return R


def get_discrete_force_Guo_without_U_experimental(i):
    """
    'Discrete lattice effects on the forcing term in the lattice Boltzmann method',  Guo et al., 2001
    version for 'Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios' A. Fakhari et. al., 2017
    """
    # first order terms only
    cs2 = 1. / 3.
    # # cs2 = Symbol('cs2')
    # eF = (ex[i] - ux) * Fx + (ey[i] - uy) * Fy  # TODO why (ey[i] - uy)
    eF = ex[i] * Fx + ey[i] * Fy
    R = w[i] * eF / (rho * cs2)
    return R


def get_discrete_force_Guo_experimental(i):
    """
    'Discrete lattice effects on the forcing term in the lattice Boltzmann method',  Guo et al., 2001
    version for 'Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios' A. Fakhari et. al., 2017
    """
    # first order terms only
    cs2 = 1. / 3.
    # # cs2 = Symbol('cs2')
    eF = (ex[i] - ux) * Fx + (ey[i] - uy) * Fy  # TODO why (ey[i] - uy)
    # eF = ex[i]*Fx + ey[i]*Fy
    R = w[i] * eF / (rho * cs2)
    return R


def get_discrete_force_Guo(i):
    """
    'Discrete lattice effects on the forcing term in the lattice Boltzmann method',  Guo et al., 2001
    version for 'Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios' A. Fakhari et. al., 2017
    """
    # extended version with second order terms
    cs2 = 1. / 3.
    temp_x = ex[i] - ux + (ex[i] * ux + ey[i] * uy) * ex[i] / cs2
    temp_y = ey[i] - uy + (ex[i] * ux + ey[i] * uy) * ey[i] / cs2
    R = w[i] * (temp_x * Fx + temp_y * Fy) / (rho * cs2)
    return R


def get_discrete_force_interface_tracking(i):
    """
    'Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios' A. Fakhari et. al., 2017
    eq7 in cm
    """
    # R = F_phi_coeff * w[i]*(ex[i]*phi_norm_grad_x + ey[i]*phi_norm_grad_y)  #improve results a bit

    # try Guo:
    # extended version with second order terms
    cs2 = 1. / 3.
    temp_x = ex[i] - ux + (ex[i] * ux + ey[i] * uy) * ex[i] / cs2
    temp_y = ey[i] - uy + (ex[i] * ux + ey[i] * uy) * ey[i] / cs2
    R = w[i] * (temp_x * F_phi_x + temp_y * F_phi_y) / cs2
    return R


def get_discrete_m(m, n, fun):
    k = 0
    for i in range(9):
        # pop = get_pop_eq(i)
        # pop = p_star * get_gamma(i)
        # pop = Symbol('f[%d]' % i)
        pop = fun(i)
        k += pow(ex[i], m) * pow(ey[i], n) * pop

    return round_and_simplify(k)


def get_discrete_cm(m, n, fun):
    k = 0
    for i in range(9):
        # pop = get_pop_eq(i)
        # pop = p_star * get_gamma(i)
        # pop = Symbol('f[%d]' % i)
        pop = fun(i)
        k += pow((ex[i] - ux), m) * pow((ey[i] - uy), n) * pop

    return round_and_simplify(k)


def get_mom_vector_from_discrete_def(fun, discrete_transform):
    #  for example: discrete_transform=get_discrete_cm

    mom_ = [discrete_transform(0, 0, fun),
            discrete_transform(1, 0, fun),
            discrete_transform(0, 1, fun),
            discrete_transform(2, 0, fun),
            discrete_transform(0, 2, fun),
            discrete_transform(1, 1, fun),
            discrete_transform(2, 1, fun),
            discrete_transform(1, 2, fun),
            discrete_transform(2, 2, fun)
           ]
    return Matrix([mom_]).transpose()


def get_mom_vector_from_shift_Mat(fun, Mat):
    pop = Matrix([fun(i) for i in range(9)])
    # pop = Matrix(9, 1, lambda i,j: i+j)  # column vect
    cm_ = Mat * pop  #for example: Mat=Nraw * Mraw)
    cm_ = round_and_simplify(cm_)
    return Matrix([cm_])


