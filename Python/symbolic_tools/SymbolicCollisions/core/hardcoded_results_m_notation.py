
from sympy.matrices import Matrix
from SymbolicCollisions.core.cm_symbols import m00, rho, cp, cs2_thermal, \
    Fx, Fy, Fz, F_phi_x, F_phi_y, F_phi_z, \
    ux, uy, ux2, uy2, uxuy, \
    uz, uz2, uxuz, uyuz

from SymbolicCollisions.core.cm_symbols import Temperature as T
from SymbolicCollisions.core.cm_symbols import Enthalpy as H
from SymbolicCollisions.core.cm_symbols import cht_gamma, Sigma2asSymbol

# save time and hardcode some of the results

hardcoded_cm_eq_cht_D3Q27 = Matrix([
    H,     # T*cp*rho,
    0,
    Sigma2asSymbol * H,     # 1/3.*T*cht_gamma,
    0,
    0,
    0,
    Sigma2asSymbol * H,     # 1/3.*T*cht_gamma,
    0,
    Sigma2asSymbol * Sigma2asSymbol * H,     # 1/3.*T*cht_gamma,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    Sigma2asSymbol * H,     # 1/9.*T*cht_gamma*cht_gamma/(cp*rho),
    0,
    Sigma2asSymbol * Sigma2asSymbol * H,     # 1/9.*T*cht_gamma*cht_gamma/(cp*rho),
    0,
    0,
    0,
    Sigma2asSymbol * Sigma2asSymbol * H,     # 1/9.*T*cht_gamma*cht_gamma/(cp*rho),
    0,
    Sigma2asSymbol * Sigma2asSymbol * Sigma2asSymbol * H,     # T*(1/27.*cht_gamma*cht_gamma*cht_gamma)/(cp*cp*rho*rho),
])
