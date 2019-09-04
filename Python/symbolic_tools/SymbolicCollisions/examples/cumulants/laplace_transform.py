from sympy import Function, Symbol, symbols, Derivative, preview, simplify, collect,  Poly
from sympy import diff, ln, pprint
import sympy as sym
import sympy as sp
from sympy import integrate, oo

from sympy.integrals.transforms import mellin_transform, laplace_transform
from sympy.abc import x, s, t
from sympy.matrices import Matrix, diag
from SymbolicCollisions.core.printers import round_and_simplify
from SymbolicCollisions.core.cm_symbols import moments_dict
from SymbolicCollisions.core.cm_symbols import s3D, s2D, s_x, u3D, u2D, dzeta3D, dzeta2D, rho, sigma2
import re

# print("---------- Normal Distribution Function and its Laplace tranform ----------")
# normal_df = 1/sp.sqrt(2*sp.pi*sigma2)
# normal_df *= sp.exp(- (x-mu)*(x-mu) / (2 * sigma2))
# transformed_normal_df = laplace_transform(normal_df, x, s)
# transformed_normal_df = round_and_simplify(transformed_normal_df)
# print(transformed_normal_df)
# # https://math.stackexchange.com/questions/1540880/laplace-transform-of-normal-distribution-function
# #(-0.353553390593274*sqrt(2)*rho*(erf(sqrt(2)*(Sigma*s - u)/(2*sqrt(Sigma))) - 1)*exp(s*(Sigma*s/2 - u)), -oo, True)
#
# print("---------- Example: Laplace transform of cos(a*t) ----------")
# a = Symbol('a', positive=True)
# fun = sp.cos(a*t)
# Lf = laplace_transform(fun, t, s)
# print(Lf)



def get_Maxwellian_DF_v2(psi, u, sigma2, dzeta):
    """
    :param u: velocity (x,y,z)
    :param sigma2: variance of the distribution
    :param psi: quantity of interest aka scaling function like density
    :return: continuous, local Maxwell-Boltzmann distribution
    'Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments'
    Kannan N. Premnath, Sanjoy Banerjee, 2009
    eq 22
    """

    PI = sp.pi
    dzeta_minus_u = dzeta - u
    dzeta_u2 = dzeta_minus_u.dot(dzeta_minus_u)

    # thank you sympy...
    # hacks:
    # for 2D
    # df = psi / pow(2 * PI * self.cs2, 2/2)
    # LOL: 2/2 gives not simplified result for m22 on d2q9:
    # 1.0 * m00 * (RT * u.y ** 2 - RT ** 1.0 * u.y ** 2 + RT ** 2.0);
    # while 1 does ;p
    # RT ** 2 * m00;

    dim = len(dzeta)  # number od dimensions
    # df = psi / pow(2 * PI * self.cs2, dim / 2)  # this is to difficult for sympy :/
    if dim == 1:
        df = psi / sp.sqrt(2 * PI * sigma2)
    elif dim == 2:
        df = psi / (2 * PI * sigma2)  # 2D version hack
    elif dim == 3:
        df = psi / pow(2 * PI * sigma2, dim / 2)  # this may be to difficult for sympy :/
        # if self.cs2 != 1. / 3.:
        #     warnings.warn("Sympy may have problem with 3D non isothermal version (cs2=RT) \n "
        #                   "It also can't simplify it, thus check the raw output", UserWarning)
    else:
        raise Exception(f"Wrong dimension: {dim}")

    df *= sp.exp(-dzeta_u2 / (2 * sigma2))
    return df


def my_laplace_transform_of_edf(_fun, _dzeta, _u, _s):
    # _fun = test_get_Maxwellian_DF_v2(rho, _u, sigma2, _dzeta)
    _fun *= sp.exp(-_s.dot(_dzeta))
    lim = [(dim, -oo, oo) for dim in _dzeta]
    # lim = [(dim, 0, oo) for dim in _dzeta]  # TODO: sympy's transform is from 0 to inf
    # lim = [(dzeta_x, -oo, oo), (dzeta_y, -oo, oo), (dzeta_z, -oo, oo)]  # oj, sympy w 3D tym razem nie scalkuje ;p
    result = integrate(_fun, *lim)
    return round_and_simplify(result)


def calc_cumulant(_fun, direction, order):
    cgf = ln(_fun)
    deriv = Derivative(cgf, (direction, order))
    all_terms = deriv.doit()
    all_terms = simplify(all_terms)
    return all_terms


# print("---------- this works in 1D... ----------")
# mb_edf_1D = get_Maxwellian_DF_v2(rho, u1D, sigma2, dzeta1D)
# l_eqdf = laplace_transform(mb_edf_1D, dzeta_x, s_x)
# print(f"sympy's 1D l_edf = \n{l_eqdf}")

print("---------- it shall work in 3D as well... ----------")
mb_edf = get_Maxwellian_DF_v2(rho, u2D, sigma2, dzeta2D)
my_l_edf = my_laplace_transform_of_edf(mb_edf, dzeta2D, u2D, s2D)
print(f"my_l_edf = {my_l_edf}")
c = calc_cumulant(my_l_edf, order=0, direction=s_x)
print(f"derivative is \n {c}")  # the the cumulant is at s = 0
result = c.subs(s_x, 0)
print(f"substituting s = 0 to get cumulant:\n {result}")  # the the cumulant is at s = 0


cgf = ln(my_l_edf)
md = moments_dict['D2Q9']

for mno in md:
    d = cgf
    for s_i, mno_i in zip(s3D, mno):
        d = Derivative(d, (s_i, mno_i), evaluate=True)
    d_at_0 = simplify(d)

    for s_i in s3D:
        d_at_0 = d_at_0.subs(s_i, 0)

    s_mno = re.sub(r',', '', str(mno))
    s_mno = re.sub(r' ', '', s_mno)
    s_mno = re.sub(r'\(', '', s_mno)
    s_mno = re.sub(r'\)', '', s_mno)

    print(f"c{s_mno} =  {d_at_0}")  # the the cumulant is at s = 0


print("try another script")

from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import F2D, F3D, e_D2Q9, e_D3Q27
from SymbolicCollisions.core.printers import print_as_vector
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho, sigma2)
lattice = 'D2Q9'
cum_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cumulants,
                                           moments_order=moments_dict[lattice],
                                           serial_run=False)
print_as_vector(cum_eq, outprint_symbol='c', output_order_of_moments=moments_dict[lattice])
print("bye")
