
from sympy.matrices import Matrix
from SymbolicCollisions.core.cm_symbols import m00, rho, cs2_thermal, \
    Fx, Fy, Fz, F_phi_x, F_phi_y, F_phi_z, \
    ux, uy, ux2, uy2, uxuy


# save time and hardcode some of the results
hardcoded_F_cm_hydro_density_based_D2Q9 = Matrix([
    0,
    Fx * m00 / rho,
    Fy * m00 / rho,
    0,
    0,
    0,
    Fy * m00 / (rho * 3.),
    Fx * m00 / (rho * 3.),
    0,
])

hardcoded_F_cm_hydro_density_based_D3Q19 = Matrix([
    0,
    Fx*m00/rho,
    Fy*m00/rho,
    Fz*m00/rho,
    0,
    0,
    0,
    0,
    0,
    0,
    1/3.*Fx*m00/rho,
    1/3.*Fx*m00/rho,
    1/3.*Fy*m00/rho,
    1/3.*Fz*m00/rho,
    1/3.*Fy*m00/rho,
    1/3.*Fz*m00/rho,
    0,
    0,
    0,
])

hardcoded_F_cm_pf_D3Q19 = Matrix([
    0,
    F_phi_x,
    F_phi_y,
    F_phi_z,
    0,
    0,
    0,
    0,
    0,
    0,
    1/3.*F_phi_x,
    1/3.*F_phi_x,
    1/3.*F_phi_y,
    1/3.*F_phi_z,
    1/3.*F_phi_y,
    1/3.*F_phi_z,
    0,
    0,
    0,
])


hardcoded_F_cm_He_hydro_LB_velocity_based_D2Q9 = Matrix([
    0,
    Fx / rho,
    Fy / rho,
    0,
    0,
    0,
    Fy / (rho * 3.),
    Fx / (rho * 3.),
    0,
])

hardcoded_F_cm_Guo_hydro_LB_velocity_based_D2Q9 = Matrix([
    0,
    Fx / rho,
    Fy / rho,
    0,
    0,
    0,
    (-2.0*Fx*uxuy - Fy*ux2 + 1./3.*Fy) / rho,
    (-Fx*uy2 + 1./3.*Fx - 2.0*Fy*uxuy) / rho,
    4.0 * uxuy * (Fx * uy + Fy * ux) / rho,
])

hardcoded_F_cm_pf_D2Q9 = Matrix([
    0,
    F_phi_x,
    F_phi_y,
    0,
    0,
    0,
    F_phi_y / 3.,
    F_phi_x / 3.,
    0,
])


hardcoded_cm_eq_compressible_D2Q9 = Matrix([
    m00,
    0,
    0,
    1 / 3. * m00,
    1 / 3. * m00,
    0,
    0,
    0,
    1. / 9. * m00,
    ])

hardcoded_cm_eq_compressible_D2Q9_thermal = Matrix([
    m00,
    0,
    0,
    cs2_thermal * m00,
    cs2_thermal * m00,
    0,
    0,
    0,
    cs2_thermal**2 * m00,
    ])

# order of 3D (central) moments as in
# `Three-dimensional cascaded lattice Boltzmann method:
# Improved implementation and consistent forcing scheme`
# by Linlin Fei, Kai H.  Luo,  Qing Li. 2018
hardcoded_cm_eq_compressible_D3Q19 = Matrix([
    m00,
    0,
    0,
    0,
    0,
    0,
    0,
    1 / 3. * m00,
    1 / 3. * m00,
    1 / 3. * m00,
    0,
    0,
    0,
    0,
    0,
    0,
    1 / 9. * m00,
    1 / 9. * m00,
    1 / 9. * m00,
])

hardcoded_cm_eq_incompressible_D2Q9 = Matrix([
    m00,
    -ux*(m00 - 1),
    uy*(-m00 + 1),
    m00*ux2 + 1/3.*m00 - ux2,
    m00*uy2 + 1/3.*m00 - uy2,
    uxuy*(m00 - 1.0),
    uy*(-m00 * ux2 - 1/3.*m00 + ux2 + 1/3.),
    ux*(-m00 * uy2 - 1/3.*m00 + uy2 + 1/3.),
    m00*ux2*uy2 + 1/3.*m00*ux2 + 1/3.*m00*uy2 + 1/9.*m00 - ux2*uy2 - 1/3.*ux2 - 1/3.*uy2,
    ])

# cm_eq[0] = m00;
# cm_eq[1] = u.x*(-m00 + 1);
# cm_eq[2] = u.y*(-m00 + 1);
# cm_eq[3] = m00*ux2 + 1./3.*m00 - ux2;
# cm_eq[4] = m00*uy2 + 1./3.*m00 - uy2;
# cm_eq[5] = uxuy*(m00 - 1.0);
# cm_eq[6] = u.y*(-m00*ux2 - 1./3.*m00 + ux2 + 1./3.);
# cm_eq[7] = u.x*(-m00*uy2 - 1./3.*m00 + uy2 + 1./3.);
# cm_eq[8] = m00*ux2*uy2 + 1./3.*m00*ux2 + 1./3.*m00*uy2 + 1./9.*m00 - ux2*uy2 - 1./3.*ux2 - 1./3.*uy2;
