from sympy import Symbol
from sympy.interactive.printing import init_printing
from sympy.matrices import Matrix, diag

# from sympy.matrices import eye, zeros, ones, GramSchmidt
# from sympy import simplify, Float, preorder_traversal
# import numpy as np

from importlib import import_module

init_printing(use_unicode=False, wrap_line=False, no_global=True)


def dynamic_import(abs_module_path, class_name):
    module_object = import_module(abs_module_path)
    target_class = getattr(module_object, class_name)
    return target_class


####################################################### SYMBOLS #######################################################
cs2_thermal = Symbol('RT', positive=True)  # variance of the distribution
cs2 = 1/3.  # variance of the distribution
sigma2 = Symbol('Sigma2', positive=True)  # variance of the distribution

ux = Symbol('u.x')  # don't set real=True for velocity as it freezes the test suite :/
uy = Symbol('u.y')
uz = Symbol('u.z')
u1D = Matrix([ux])
u2D = Matrix([ux, uy])
u3D = Matrix([ux, uy, uz])

ux2 = Symbol('ux2')
uy2 = Symbol('uy2')
uz2 = Symbol('uz2')

ux3 = Symbol('ux3')
uy3 = Symbol('uy3')
uz3 = Symbol('uy3')

uxuy3 = Symbol('uxuy3')

uxuy = Symbol('uxuy')
uxuz = Symbol('uxuz')
uyuz = Symbol('uyuz')

m00 = Symbol('m00', positive=True)
rho = Symbol('rho', positive=True)
Temperature = Symbol('T', positive=True)
Enthalpy = Symbol('H', positive=True)  # consider it as enthalpy
cp = Symbol('cp', positive=True)
cht_gamma = Symbol('h_stability_enhancement', positive=True)  # magic stability enhancement
# (h_stability_enhancement * 1. / 3.) / (cp * rho);

Sigma2 = cs2*cht_gamma/(cp*rho)
Sigma2asSymbol = Symbol('Sigma2', positive=True)

dzeta_x = Symbol('dzeta_x', real=True)
dzeta_y = Symbol('dzeta_y', real=True)
dzeta_z = Symbol('dzeta_z', real=True)
dzeta1D = Matrix([dzeta_x])
dzeta2D = Matrix([dzeta_x, dzeta_y])
dzeta3D = Matrix([dzeta_x, dzeta_y, dzeta_z])

s_x = Symbol('s_x', real=False)
s_y = Symbol('s_y', real=False)
s_z = Symbol('s_z', real=False)
s1D = Matrix([s_x])
s2D = Matrix([s_x, s_y])
s3D = Matrix([s_x, s_y, s_z])

Force_str = "F"
Fx = Symbol(Force_str + '.x')
Fy = Symbol(Force_str + '.y')
Fz = Symbol(Force_str + '.z')

F2D = Matrix([Fx, Fy])
F3D = Matrix([Fx, Fy, Fz])

F_phi_x = Symbol(Force_str + '.x')
F_phi_y = Symbol(Force_str + '.y')
F_phi_z = Symbol(Force_str + '.z')

omega_ade = Symbol('omega_ade')
omega_v = Symbol('omega_nu')
omega_b = Symbol('omega_bulk')  # omega_bulk='1.0/(3*bulk_visc+0.5)'

####################################################### END OF SYMBOLS #######################################################
# D2Q9 notation from TCLB
ex_D2Q9 = Matrix([0, 1, 0, -1, 0, 1, -1, -1, 1])
ey_D2Q9 = Matrix([0, 0, 1, 0, -1, 1, 1, -1, -1])
ez_D2Q9 = Matrix([0, 0, 0, 0, 0, 0, 0, 0, 0])
e_D2Q9 = ex_D2Q9.col_insert(1, ey_D2Q9)

# D3Q7 notation from TCLB
ex_D3Q7 = Matrix([0, 1, -1, 0, 0, 0, 0])
ey_D3Q7 = Matrix([0, 0, 0, 1, -1, 0, 0])
ez_D3Q7 = Matrix([0, 0, 0, 0, 0, 1, -1])

e_D3Q7 = ex_D3Q7.col_insert(1, ey_D3Q7)
e_D3Q7 = e_D3Q7.col_insert(2, ez_D3Q7)

S_relax_ADE_D3Q7 = diag(1, omega_ade, omega_ade, omega_ade, 1, 1, 1)

# D3Q15 - notation from 'LBM Principles and Practise' Book p. 89
# ex_D3Q15 = Matrix([0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1])
# ey_D3Q15 = Matrix([0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1])
# ez_D3Q15 = Matrix([0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1])

# D3Q15 - notation for phase-field as in TCLB's d3q27_pf_velocity model
ex_D3Q15 = Matrix([0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1])
ey_D3Q15 = Matrix([0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1])
ez_D3Q15 = Matrix([0, 0, 0, 0, 0, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1])

S_relax_ADE_D3Q15 = diag(1, omega_ade, omega_ade, omega_ade, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)


# D3Q19 -
# as in TCLB or '3D CLBM: Improved implementation and consistent forcing scheme' by L. Fei and Q. li 2018
# (differs from 'LBM Principles and Practise' Book)
ex_D3Q19 = Matrix([0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0])
ey_D3Q19 = Matrix([0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1])
ez_D3Q19 = Matrix([0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1])

S_relax_ADE_D3Q19 = diag(1, omega_ade, omega_ade, omega_ade, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# FROM TCLB
# d3q19 = matrix(c(
#          0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,
#          0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,
#          0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1), 19, 3);
#
# attr(d3q19,"MAT") = matrix(c(
#          1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
#        -30,-11,-11,-11,-11,-11,-11,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,
#         12, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
#          0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,
#          0, -4,  4,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,
#          0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,
#          0,  0,  0, -4,  4,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,
#          0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,
#          0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,
#          0,  2,  2, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2,
#          0, -4, -4,  2,  2,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2,
#          0,  0,  0,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,
#          0,  0,  0, -2, -2,  2,  2,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,
#          0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,
#          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,
#          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,
#          0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1, -1,  1, -1,  1,  0,  0,  0,  0,
#          0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  0,  0,  0,  0,  1, -1,  1, -1,
#          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1,  1), 19, 19);



# D3Q27 - notation for hydrodynamics as in TCLB's d3q27_pf_velocity model
ex_D3Q27 = Matrix([0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0])
ey_D3Q27 = Matrix([0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1])
ez_D3Q27 = Matrix([0, 0, 0, 0, 0, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1])
e_D3Q27 = ex_D3Q27.col_insert(1, ey_D3Q27)
e_D3Q27 = e_D3Q27.col_insert(2, ez_D3Q27)

# D3Q27 - notation from "Three-dimensional cascaded lattice
# % Boltzmann method: improved implementation and consistent forcing scheme" by Linlin Fei
# ex=[0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1];
# ey=[0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1,-1];
# ez=[0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1];


# phi_norm_grad_x = Symbol('norm_grad_phi.x')  # normalized gradient of the phase field in the x direction
# phi_norm_grad_y = Symbol('norm_grad_phi.y')  # normalized gradient of the phase field in the y direction
# phi_norm_grad_z = Symbol('norm_grad_phi.z')  # normalized gradient of the phase field in the y direction
#
# F_phi_coeff = Symbol('F_phi_coeff')  # F_phi_coeff=(1.0 - 4.0*(myPhaseF - pfavg)*(myPhaseF - pfavg))/inteface_width;



w_D2Q9 = Matrix([4. / 9, 1. / 9, 1. / 9, 1. / 9, 1. / 9, 1. / 36, 1. / 36, 1. / 36, 1. / 36])




# this matrix will produce raw moments (m=M*f) in the following order:
# [m00, m10, m01, m20, m02, m11, m21, m12, m22]
# "Modelling incompressible thermal flows using a central-moments-based lattice Boltzmann method" L. Fei et al. 2017
Mraw_D2Q9 = Matrix([
    [1, 1, 1,  1,  1, 1,  1,  1,  1],
    [0, 1, 0, -1,  0, 1, -1, -1,  1],
    [0, 0, 1,  0, -1, 1,  1, -1, -1],
    [0, 1, 0,  1,  0, 1,  1,  1,  1],
    [0, 0, 1,  0,  1, 1,  1,  1,  1],
    [0, 0, 0,  0,  0, 1, -1,  1, -1],
    [0, 0, 0,  0,  0, 1,  1, -1, -1],
    [0, 0, 0,  0,  0, 1, -1, -1,  1],
    [0, 0, 0,  0,  0, 1,  1,  1,  1]
])

# eq 10.30 from The Lattice Boltzmann Method: Principles and Practice
# T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
M_ortho_GS = Matrix([
    [ 1,  1,  1,  1,  1, 1,  1,  1,  1],
    [-4, -1, -1, -1, -1, 2,  2,  2,  2],
    [ 4, -2, -2, -2, -2, 1,  1,  1,  1],
    [ 0,  1,  0, -1,  0, 1, -1, -1,  1],
    [ 0, -2,  0,  2,  0, 1, -1, -1,  1],
    [ 0,  0,  1,  0, -1, 1,  1, -1, -1],
    [ 0,  0, -2,  0,  2, 1,  1, -1, -1],
    [ 0,  1, -1,  1, -1, 0,  0,  0,  0],
    [ 0,  0,  0,  0,  0, 1, -1,  1, -1]
])

ex_Straka_d2_q5 = Matrix([0, -1, 0, 1, 0])
ey_Straka_d2_q5 = Matrix([0, 0, -1, 0, 1])

K_ortho_Straka_d2q5 = Matrix([  # in Geier's lattice numbering CSYS
    [1,  0,  0,  4,  0],  # 0
    [1, -1,  0, -1, -1],  # 1
    [1,  0, -1, -1,  1],  # 2
    [1,  1,  0, -1, -1],  # 3
    [1,  0,  1, -1,  1],  # 4
])

Shift_ortho_Straka_d2q5 = Matrix([  # in Geier's lattice numbering CSYS
    [2,  0, 0, 0],  # 1
    [0,  2, 0, 0],  # 2
    [-4 * ux,  0,  -2, -2],  # 3
    [0,  -4 * uy,  -2,  2],  # 4
])

ex_Geier = Matrix([0, -1, -1, -1, 0, 1, 1, 1, 0])
ey_Geier = Matrix([0, 1, 0, -1, -1, -1, 0, 1, 1])

K_ortho_Geier = Matrix([  # in Geier's lattice numbering CSYS
    [1,  0,  0, -4,  0,  0,  0,  0,  4],  # 0
    [1, -1,  1,  2,  0,  1, -1,  1,  1],  # 1
    [1, -1,  0, -1,  1,  0,  0, -2, -2],  # 2
    [1, -1, -1,  2,  0, -1,  1,  1,  1],  # 3
    [1,  0, -1, -1, -1,  0, -2,  0, -2],  # 4
    [1,  1, -1,  2,  0,  1,  1, -1,  1],  # 5
    [1,  1,  0, -1,  1,  0,  0,  2, -2],  # 6
    [1,  1,  1,  2,  0, -1, -1, -1,  1],  # 7
    [1,  0,  1, -1, -1,  0,  2,  0, -2],  # 8
])

Shift_ortho_Geier = Matrix([
    [                  6,               2,             0,      0,      0, 0],  # noqa
    [                  6,              -2,             0,      0,      0, 0],  # noqa
    [                  0,               0,            -4,      0,      0, 0],  # noqa
    [            -6 * uy,         -2 * uy,        8 * ux,     -4,      0, 0],  # noqa
    [            -6 * ux,          2 * ux,        8 * uy,      0,     -4, 0],  # noqa
    [8 + 6 * (ux2 + uy2), 2 * (uy2 - ux2), -16 * ux * uy, 8 * uy, 8 * ux, 4],
])

# raw moments - interpretation
# real_t m00 = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7]  + f[8];  // m00 - m0: density
# real_t m10 =        f[1]        - f[3]        + f[5] - f[6] - f[7]  + f[8];  // m10 - m1: x momentum flux
# real_t m01 =             + f[2]        - f[4] + f[5] + f[6] - f[7]  - f[8];  // m01 - m2: y momentum flux
# real_t m20 =        f[1]        + f[3]        + f[5] + f[6] + f[7]  + f[8];  // m20 - m3
# real_t m02 =               f[2]        + f[4] + f[5] + f[6] + f[7]  + f[8];  // m02 - m4
# real_t m11 =                                    f[5] - f[6] + f[7]  - f[8];  // m11 - m5: stress tensor xy (off-diagonal)
# real_t m21 =                                    f[5] + f[6] - f[7]  - f[8];  // m21 - m6
# real_t m12 =                                    f[5] + f[6] - f[7]  - f[8];  // m12 - m7
# real_t m22 =                                    f[5] + f[6] + f[7]  + f[8];  // m22 - m8

# SHIFT MATRIX
# "Modelling incompressible thermal flows using a central-moments-based lattice Boltzmann method" L. Fei et al. 2017
NrawD2Q9 = Matrix([
    [                1,                 0,                 0,       0,       0,            0,       0,       0, 0],
    [              -ux,                 1,                 0,       0,       0,            0,       0,       0, 0],
    [              -uy,                 0,                 1,       0,       0,            0,       0,       0, 0],
    [          ux * ux,           -2 * ux,                 0,       1,       0,            0,       0,       0, 0],
    [          uy * uy,                 0,           -2 * uy,       0,       1,            0,       0,       0, 0],
    [          ux * uy,               -uy,               -ux,       0,       0,            1,       0,       0, 0],
    [    -ux * ux * uy,       2 * ux * uy,           ux * ux,     -uy,       0,      -2 * ux,       1,       0, 0],
    [    -uy * uy * ux,           uy * uy,       2 * ux * uy,       0,     -ux,      -2 * uy,       0,       1, 0],
    [ux * ux * uy * uy, -2 * ux * uy * uy, -2 * uy * ux * ux, uy * uy, ux * ux,  4 * ux * uy, -2 * uy, -2 * ux, 1],
])

# RELAXATION MATRIX
s_plus_D2Q9 = (omega_b + omega_v) / 2
s_minus_D2Q9 = (omega_b - omega_v) / 2

S_relax_hydro_D2Q9 = diag(1, 1, 1, s_plus_D2Q9, s_plus_D2Q9, omega_v, 1, 1, 1)
S_relax_hydro_D2Q9[3, 4] = s_minus_D2Q9
S_relax_hydro_D2Q9[4, 3] = s_minus_D2Q9

S_relax_ADE_D2Q9 = diag(1, omega_ade, omega_ade, 1, 1, 1, 1, 1, 1)

# Both S_relax_D3Q27 and order of 3D (central) moments as in
# `Three-dimensional cascaded lattice Boltzmann method:
# Improved implementation and consistent forcing scheme`
# by Linlin Fei, Kai H.  Luo,  Qing Li. 2018

# bulk_visc_3D = cs2*(2/3)*(1/sb -0.5)  # this is different than in 2D for some reason...
#  kin_visc_3D = cs2(1/s_v - 0.5)

s_plus_D3Q27 = (omega_b + 2*omega_v)/3
s_minus_D3Q27 = (omega_b - omega_v)/3
S_relax_hydro_D3Q27 = diag(1, 1, 1, 1, omega_v, omega_v, omega_v, s_plus_D3Q27, s_plus_D3Q27, s_plus_D3Q27, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
S_relax_hydro_D3Q27[7, 8] = s_minus_D3Q27
S_relax_hydro_D3Q27[7, 9] = s_minus_D3Q27
S_relax_hydro_D3Q27[8, 7] = s_minus_D3Q27
S_relax_hydro_D3Q27[8, 9] = s_minus_D3Q27
S_relax_hydro_D3Q27[9, 7] = s_minus_D3Q27
S_relax_hydro_D3Q27[9, 8] = s_minus_D3Q27

# S_relax_ADE_D3Q27 = diag(1, omega_ade, omega_ade, omega_ade, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
S_relax_ADE_D3Q27 = diag(1, omega_ade, omega_ade, omega_ade, 1, 1, 1, 1, 1, 1, omega_ade, omega_ade, omega_ade, omega_ade, omega_ade, omega_ade, 1, 1, 1, 1, 1, 1, 1, omega_ade, omega_ade, omega_ade, 1)

S_relax_MRT_GS = diag(1, 1, 1, 1, 1, 1, 1, omega_v, omega_v)  #
# S_relax_MRT_GS = diag(0, 0, 0, 0, 0, 0, 0, sv, sv)   #


moments_dict = {
    # order of 2D (central) moments as in
    # `Modeling incompressible thermal flows using a central-moments-based lattice Boltzmann method`
    # by Linlin Fei, Kai Hong Luo Chuandong Lin , Qing Li. 2017
    'D2Q5': [(0, 0, 0),
             (1, 0, 0),
             (0, 1, 0),
             (2, 0, 0),
             (0, 2, 0)],

    'D2Q9': [(0, 0, 0),
             (1, 0, 0),
             (0, 1, 0),
             (2, 0, 0),
             (0, 2, 0),
             (1, 1, 0),
             (2, 1, 0),
             (1, 2, 0),
             (2, 2, 0)
             ],

    # order of 3D (central) moments as in
    # `Three-dimensional cascaded lattice Boltzmann method:
    # Improved implementation and consistent forcing scheme`
    # by Linlin Fei, Kai H.  Luo,  Qing Li. 2018
    'D3Q7': [(0, 0, 0),
             (1, 0, 0),
             (0, 1, 0),
             (0, 0, 1),
             (2, 0, 0),
             (0, 2, 0),
             (0, 0, 2)],

    'D3Q15': [(0, 0, 0),
              (1, 0, 0),
              (0, 1, 0),
              (0, 0, 1),
              (2, 0, 0),
              (0, 2, 0),
              (0, 0, 2),
              (1, 1, 1),  # skipped in      D3Q19, D3Q7
              (2, 1, 1),  # skipped in      D3Q19, D3Q7
              (1, 2, 1),  # skipped in      D3Q19, D3Q7
              (1, 1, 2),  # skipped in      D3Q19, D3Q7
              (1, 2, 2),  # skipped in      D3Q19, D3Q7
              (2, 1, 2),  # skipped in      D3Q19, D3Q7
              (2, 2, 1),  # skipped in      D3Q19, D3Q7
              (2, 2, 2),  # skipped in      D3Q19, D3Q7
              ],

    'D3Q19': [(0, 0, 0),
              (1, 0, 0),
              (0, 1, 0),
              (0, 0, 1),
              (1, 1, 0),  # skipped in D3Q15, D3Q7
              (1, 0, 1),  # skipped in D3Q15, D3Q7
              (0, 1, 1),  # skipped in D3Q15, D3Q7
              (2, 0, 0),
              (0, 2, 0),
              (0, 0, 2),
              (1, 2, 0),  # skipped in D3Q15, D3Q7
              (1, 0, 2),  # skipped in D3Q15, D3Q7
              (2, 1, 0),  # skipped in D3Q15, D3Q7
              (2, 0, 1),  # skipped in D3Q15, D3Q7
              (0, 1, 2),  # skipped in D3Q15, D3Q7
              (0, 2, 1),  # skipped in D3Q15, D3Q7
              (2, 2, 0),  # skipped in D3Q15, D3Q7
              (2, 0, 2),  # skipped in D3Q15, D3Q7
              (0, 2, 2),  # skipped in D3Q15, D3Q7
              ],

    'D3Q27': [(0, 0, 0),
              (1, 0, 0),
              (0, 1, 0),
              (0, 0, 1),
              (1, 1, 0),  # skipped in D3Q15, D3Q7
              (1, 0, 1),  # skipped in D3Q15, D3Q7
              (0, 1, 1),  # skipped in D3Q15, D3Q7
              (2, 0, 0),
              (0, 2, 0),
              (0, 0, 2),
              (1, 2, 0),  # skipped in D3Q15, D3Q7
              (1, 0, 2),  # skipped in D3Q15, D3Q7
              (2, 1, 0),  # skipped in D3Q15, D3Q7
              (2, 0, 1),  # skipped in D3Q15, D3Q7
              (0, 1, 2),  # skipped in D3Q15, D3Q7
              (0, 2, 1),  # skipped in D3Q15, D3Q7
              (1, 1, 1),  # skipped in      D3Q19, D3Q7
              (2, 2, 0),  # skipped in D3Q15, D3Q7
              (2, 0, 2),  # skipped in D3Q15, D3Q7
              (0, 2, 2),  # skipped in D3Q15, D3Q7
              (2, 1, 1),  # skipped in      D3Q19, D3Q7
              (1, 2, 1),  # skipped in      D3Q19, D3Q7
              (1, 1, 2),  # skipped in      D3Q19, D3Q7
              (1, 2, 2),  # skipped in      D3Q19, D3Q7
              (2, 1, 2),  # skipped in      D3Q19, D3Q7
              (2, 2, 1),  # skipped in      D3Q19, D3Q7
              (2, 2, 2),  # skipped in      D3Q19, D3Q7
              ],

    # beware that, R in tclb generates different order of moments and lattice vectors (e),
    # 'D3Q27_r': [(0, 0, 0),
    #           (1, 0, 0),
    #           (2, 0, 0),
    #           (0, 1, 0),
    #           (1, 1, 0),
    #           (2, 1, 0),
    #           (0, 2, 0),
    #           (1, 2, 0),
    #           (2, 2, 0),
    #           (0, 0, 1),
    #           (1, 0, 1),
    #           (2, 0, 1),
    #           (0, 1, 1),
    #           (1, 1, 1),
    #           (2, 1, 1),
    #           (0, 2, 1),
    #           (1, 2, 1),
    #           (2, 2, 1),
    #           (0, 0, 2),
    #           (1, 0, 2),
    #           (2, 0, 2),
    #           (0, 1, 2),
    #           (1, 1, 2),
    #           (2, 1, 2),
    #           (0, 2, 2),
    #           (1, 2, 2),
    #           (2, 2, 2),
    #           ],
}
