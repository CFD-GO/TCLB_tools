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
from SymbolicCollisions.core.cm_symbols import ex, ey, ux, uy, w, m00, \
    Fx, Fy, F_phi_x, F_phi_y, rho, dzeta_x, dzeta_y, \
    Nraw, Mraw, M_ortho_GS


from SymbolicCollisions.core.printers import round_and_simplify

init_printing(use_unicode=False, wrap_line=False, no_global=True)


def get_DF(print_symbol='default_symbol2', start=0, end=9):
    symbols_ = [Symbol("%s[%d]" % (print_symbol, i)) for i in range(start, end)]

    return Matrix(symbols_)


def get_m00(print_symbol='default_symbol3'):
    m00_ = Symbol("%s[%d]" % (print_symbol, 0))

    for i in range(1, 9):
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


def get_discrete_force_Guo_first_order(i):
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


def get_discrete_force_Guo_second_order(i):
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


def get_continuous_weight(dzeta=(dzeta_x, dzeta_y)):

    """
    PhD Thesis: `The lattice Boltzmann method: Fundamentals and acoustics`
    by Erlend Magnus Viggen
    4.1  The discrete-velocity Boltzmann equation, pp75
    :param i: i-th lattice direction
    :return: returns weight in i-th lattice direction
    """
    e2 = dzeta[0] * dzeta[0] + dzeta[1] * dzeta[1]

    cs2 = 1. / 3.
    dim = 2  # dimension of the space
    w_ = 1. / pow((2 * pi * cs2), dim / 2.)
    w_ *= exp(-e2 / (2 * cs2))
    return w_


def get_continuous_force_Guo(dzeta=(dzeta_x, dzeta_y)):
    cs2 = 1. / 3.
    # cs2 = Symbol('cs2')
    # extended version with second order terms
    temp_x = dzeta[0] - ux + (dzeta[0] * ux + dzeta[1] * uy) * dzeta[0] / cs2
    temp_y = dzeta[1] - uy + (dzeta[0] * ux + dzeta[1] * uy) * dzeta[1] / cs2

    result = get_continuous_weight(dzeta) * (temp_x * Fx + temp_y * Fy) / (rho * cs2)
    return result


def get_continuous_Maxwellian_DF(dzeta=(dzeta_x, dzeta_y), psi=m00, u=(ux, uy)):
    """
    :param dzeta: direction (x,y)
    :param u: velocity (x,y)
    :param psi: quantity of interest aka scaling function like density
    :return: continuous, local Maxwell-Boltzmann distribution
    'Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments'
    Kannan N. Premnath, Sanjoy Banerjee, 2009
    eq 22
    """

    cs2 = 1. / 3.
    dzeta_u2 = 0
    for dzeta_i, u_i in zip(dzeta, u):
        dzeta_u2 += (dzeta_i-u_i)*(dzeta_i-u_i)

    DF = psi / (2 * pi * cs2)
    DF *= exp(-dzeta_u2 / (2 * cs2))

    return DF

def get_continuous_hydro_DF(dzeta=(dzeta_x, dzeta_y)):

    DF_p = get_continuous_Maxwellian_DF(dzeta=(dzeta_x, dzeta_y),psi=(m00-1), u=(0, 0))
    DF_gamma = get_continuous_Maxwellian_DF(dzeta=(dzeta_x, dzeta_y),psi=(1), u=(ux, uy))
    return DF_p + DF_gamma


def get_continuous_force_He_hydro_DF(dzeta=(dzeta_x, dzeta_y)):
    """
    'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
    """
    cs2 = 1. / 3.
    eu = dzeta[0] * Fx + dzeta[1] * Fy
    dzeta_2 = dzeta[0]*dzeta[0] + dzeta[1] * dzeta[1]
    DF_p = (m00-1) / (2 * pi * cs2)*exp(-dzeta_2 / (2 * cs2))
    # DF_p = get_continuous_Maxwellian_DF(dzeta=dzeta, psi=(m00-1), u=(0, 0))

    euF = (dzeta[0] - ux) * Fx + (dzeta[1] - uy) * Fy
    dzeta_u_2 = (dzeta[0] - ux) * (dzeta[0] - ux) + (dzeta[1] - uy) * (dzeta[1] - uy)
    DF_gamma = 1 / (2 * pi * cs2)*exp(-dzeta_u_2 / (2 * cs2))
    # DF_gamma = get_continuous_Maxwellian_DF(dzeta=dzeta, psi=1, u=(ux, uy))

    R = -(eu * DF_p + euF * DF_gamma)/(rho * cs2)
    R = -R  # `-` sign is skipped to ease code copy-paste ;p
    return R


def get_continuous_force_He_MB(dzeta=(dzeta_x, dzeta_y)):
    """
    'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
    Use Maxwellian to calculate equilibria
    """
    cs2 = 1. / 3.
    # cs2 = Symbol('cs2')
    euF = (dzeta[0] - ux) * Fx + (dzeta[1] - uy) * Fy
    R = get_continuous_Maxwellian_DF(dzeta) * euF / (rho * cs2)
    return R


def get_continuous_m(m, n, DF):
    fun = DF((dzeta_x, dzeta_y)) * pow((dzeta_x), m) * pow((dzeta_y), n)

    result = integrate(fun, (dzeta_x, -oo, oo), (dzeta_y, -oo, oo))
    return round_and_simplify(result)


def get_continuous_cm(m, n, DF):
    fun = DF((dzeta_x, dzeta_y)) * pow((dzeta_x - ux), m) * pow((dzeta_y - uy), n)

    result = integrate(fun, (dzeta_x, -oo, oo), (dzeta_y, -oo, oo))
    return round_and_simplify(result)


def get_mom_vector_from_continuous_def(fun, continuous_transformation):
    # for example: continous_transformation=get_continuous_cm

    m_ = [continuous_transformation(0, 0, fun),
          continuous_transformation(1, 0, fun),
          continuous_transformation(0, 1, fun),
          continuous_transformation(2, 0, fun),
          continuous_transformation(0, 2, fun),
          continuous_transformation(1, 1, fun),
          continuous_transformation(2, 1, fun),
          continuous_transformation(1, 2, fun),
          continuous_transformation(2, 2, fun)
          ]
    return Matrix([m_])

