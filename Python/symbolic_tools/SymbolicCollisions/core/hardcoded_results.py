
from sympy.matrices import Matrix
from SymbolicCollisions.core.cm_symbols import m00, rho, \
    Fx, Fy, F_phi_x, F_phi_y, \
    ux, uy, ux2, uy2, uxuy


# save time and hardcode some of the results
hardcoded_F_cm_hydro_LB_density_based = Matrix([
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

hardcoded_F_cm_He_hydro_LB_velocity_based = Matrix([
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

hardcoded_F_cm_Guo_hydro_LB_velocity_based = Matrix([
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

hardcoded_F_cm_pf = Matrix([
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



hardcoded_cm_pf_eq = Matrix([m00,
                             0,
                             0,
                             1./3. * m00,
                             1./3. * m00,
                             0,
                             0,
                             0,
                             1./9. * m00,
                             ])



hardcoded_cm_hydro_eq = Matrix([
    m00,
    ux * (-m00 + 1),
    uy * (-m00 + 1),
    m00 * ux2 + 1. / 3. * m00 - ux2,
    m00 * uy2 + 1. / 3. * m00 - uy2,
    uxuy * (m00 - 1.0),
    uy * (-m00 * ux2 - 1. / 3. * m00 + ux2 + 1. / 3.),
    ux * (-m00 * uy2 - 1. / 3. * m00 + uy2 + 1. / 3.),
    m00 * ux2 * uy2 + 1. / 3. * m00 * ux2 + 1. / 3. * m00 * uy2 + 1. / 9. * m00 - ux2 * uy2 - 1. / 3. * ux2 - 1. / 3. * uy2,
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