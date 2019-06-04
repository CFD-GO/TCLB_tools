import numpy as np
from Benchmarks.HotKarman_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein, get_Nu_cylinder_by_Zukauskas_Jacob


# input in LB units
Re = 10
u = 0.01
D = 30

Pr = 1
rho = 1
cp = 1
# outputs

v = u*D/Re
tau = 3*v + 0.5
print(f"viscosity={v:0.6f} \t tau={tau:0.3f} \t omega={1/tau:0.3f}")

alpha = v/Pr  # thermal diffusivity
k = alpha*rho*cp  # thermal conductivity
tau_k = 3*k + 0.5
print(f"conductivity={k:0.6f} \t tau_k={tau_k:0.3f} \t omega_k={1/tau_k:0.3f}")


Nu = get_Nu_cylinder_by_Churchill_Bernstein(Re=Re, Pr=Pr)
# Nu = get_Nu_cylinder_by_Zukauskas_Jacob(Re=Re, Pr=Pr)
print(f"Nu={Nu:0.3f}")

ht_coeff_correl = Nu*k/D
print(f"average heat transfer coefficient from correlations={ht_coeff_correl:0.6f} [W/(m2*K)]")


# Newton's law of cooling
# To find heat flux in Paraview: Slice -> Calculator (T*U.X) --> Integrate Variables

Surface = np.pi*D * 1  # [m2]
# Surface = 150
# q_conv = 0.05  # 0.053 [W]
q_conv = 0.218  # [W]
T_surf = 1
T_inf = 0
ht_coeff_experimental= q_conv / (Surface * (T_surf - T_inf))
print(f"experimental heat transfer coefficient = {ht_coeff_experimental:0.6f} [W/(m2*K)]")

print(f"experimental/correlation heat transfer coefficient = {ht_coeff_experimental/ht_coeff_correl:0.6f} [-]")

# FYI
# # Air at ISA (15C)
# nu = 1.47 * 1E-5
# k = 0.025
# alpha = 2.04 * 1E-5
# rho = 1.225
# cp = 1006
# Pr=0.7323
#
# # Air at 100C
# nu = 2.306 * 1E-5
# k = 0.03095
# alpha = 3.243 * 1E-5
# rho = 0.9458
# cp = 1009
# Pr = 0.7111
#


