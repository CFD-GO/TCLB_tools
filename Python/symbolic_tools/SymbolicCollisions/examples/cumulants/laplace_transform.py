from sympy import Function, Symbol, symbols, Derivative, preview, simplify, collect,  Poly
from sympy import diff, ln, pprint
import sympy as sym
import sympy as sp
from sympy import integrate, oo

from sympy.integrals.transforms import mellin_transform, laplace_transform
from sympy.abc import x, s, t
from sympy.matrices import Matrix, diag
from SymbolicCollisions.core.printers import round_and_simplify

mu = Symbol('mu', real=True)
sigma2 = Symbol('Sigma2', positive=True)
#
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

rho = Symbol('rho', positive=True)
ux = Symbol('u.x', real=True)
uy = Symbol('u.y', real=True)
uz = Symbol('u.z', real=True)

dzeta_x = Symbol('dzeta_x', real=True)
dzeta_y = Symbol('dzeta_y', real=True)
dzeta_z = Symbol('dzeta_z', real=True)

u1D = Matrix([ux])
u2D = Matrix([ux, uy])
u3D = Matrix([ux, uy, uz])

dzeta1D = Matrix([dzeta_x])
dzeta2D = Matrix([dzeta_x, dzeta_y])
dzeta3D = Matrix([dzeta_x, dzeta_y, dzeta_z])


s_x = Symbol('s_x', real=True)  # TODO: or not real...
s_y = Symbol('s_y', real=True)
s_z = Symbol('s_z', real=True)

s1D = Matrix([s_x])
s2D = Matrix([s_x, s_y])
s3D = Matrix([s_x, s_y, s_z])


def test_get_Maxwellian_DF_v2(psi, u, sigma2, dzeta):
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

    if dim == 2:
        df = psi / (2 * PI * sigma2)  # 2D version hack
    else:
        df = psi / pow(2 * PI * sigma2, dim / 2)  # this may be to difficult for sympy :/
        # if self.cs2 != 1. / 3.:
        #     warnings.warn("Sympy may have problem with 3D non isothermal version (cs2=RT) \n "
        #                   "It also can't simplify it, thus check the raw output", UserWarning)

    df *= sp.exp(-dzeta_u2 / (2 * sigma2))
    return df


def get_laplace_transform_of_edf(_fun, _dzeta, _u, _s):
    # _fun = test_get_Maxwellian_DF_v2(rho, _u, sigma2, _dzeta)
    _fun *= sp.exp(-_s.dot(_dzeta))
    # lim = [(dim, -oo, oo) for dim in _dzeta]
    lim = [(dim, 0, oo) for dim in _dzeta]  # TODO: sympy's transform is from 0 to inf... and it works in 2D
    result = integrate(_fun, *lim)
    return round_and_simplify(result)


def calc_cumulant(_fun, direction, order):
    cgf = ln(_fun)
    deriv = Derivative(cgf, (direction, order))
    all_terms = deriv.doit()
    all_terms = simplify(all_terms)
    return all_terms

print("---------- this works in 1D... ----------")
mb_edf = test_get_Maxwellian_DF_v2(rho, u1D, sigma2, dzeta1D)
transformed_lb_eqdf = laplace_transform(mb_edf, dzeta_x, s_x)
print(transformed_lb_eqdf)

print("---------- this shall work in 3D as well... ----------")
L_edf = get_laplace_transform_of_edf(mb_edf, dzeta2D, u2D, s2D)
print("L_edf = " + str(L_edf))
c = calc_cumulant(L_edf, order=1, direction=s_x)
result = c.subs(s_x, 0)
print("cumulant is at s = 0:\n" + str(result))  # the the cumulant is at s = 0

print("bye")


# lb_edf = test_get_Maxwellian_DF_v2(rho, u2D, sigma2, dzeta2D)
# transformed_lb_eqdf = laplace_transform(lb_edf, (dzeta_x, dzeta_y), s )
#
# print(transformed_lb_eqdf)
#
# cgf = ln(lb_edf)
# deriv = Derivative(cgf, (dzeta, 1))
# all_terms = deriv.doit()
# all_terms = simplify(all_terms)
#
# print("cumulant is at s = 0:\n" + str(all_terms))  # the the cumulant is at s = 0