import numpy as np
from Benchmarks.HT_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein


# input in LB units
Re = 100
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
k = alpha*rho*cp # thermal conductivity
tau_k = 3*k + 0.5
print(f"conductivity={k:0.6f} \t tau_k={tau_k:0.3f} \t omega_k={1/tau_k:0.3f}")


Nu = get_Nu_cylinder_by_Churchill_Bernstein(Re=Re, Pr=Pr)
print(f"Nu={Nu:0.3f}")

ht_coeff = Nu*k/D
print(f"average heat transfer coefficient from correlations={ht_coeff:0.3f} [W/(m2*K)]")


# Newton's law of cooling
Surface = np.pi*D *1 # [m2]
q_conv = 4068/150 #0.053  # [W]
T_surf = 1
T_inf = 0
ht_coeff_N=q_conv/Surface*(T_surf-T_inf)
print(f"experimental heat transfer coefficient = {ht_coeff_N:0.3f} [W/(m2*K)]")
# some notes on channel dimension

