import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

from Benchmarks.PoiseuilleFlow.PoiseuilleAnal import TwoPhasePoiseuilleAnal, calc_gx_between_plates
from Benchmarks.PoiseuilleFlow.TwoPhasePoiseuilleFD_between_plates import TwoPhasePoiseuilleFD_between_plates, get_tanh_profile

# magic_ratio = 1
# rho_h = magic_ratio * 54.6
# rho_l = 1
# kin_visc_h = 1
# kin_visc_l = magic_ratio

# rho_h = 1000
# rho_l = 1.22
# kin_visc_h = 1.0E-06
# kin_visc_l = 1.5E-05

rho_h = 1
rho_l = 1
kin_visc_h = 100
kin_visc_l = 1

mu_h = rho_h * kin_visc_h
mu_l = rho_l * kin_visc_l
mu_ratio = mu_l / mu_h

h = 500
y_ = np.linspace(-h, h, 10000)

uc = 0.01


# gx = 1
gx = calc_gx_between_plates(uc, mu_l, mu_h, rho_l, rho_h, h)
pa = TwoPhasePoiseuilleAnal(gx=gx, mu_l=mu_l, mu_h=mu_h, rho_h=rho_h, rho_l=rho_l, h=h)
u_a = np.array([pa.get_u_profile(y_[i]) for i in range(len(y_))])

p_fd = TwoPhasePoiseuilleFD_between_plates(gx=gx, mu_l=mu_l, mu_h=mu_h, rho_h=rho_h, rho_l=rho_l, r=h)
u_fd = p_fd.get_u_profile(y_, W=5)

gx_air_water = calc_gx_between_plates(uc, mu_l=1.22 * 1.5E-05, mu_h=1000 * 1.0E-06, rho_h=1000, rho_l=1.22, h=h)
pa_air_water = TwoPhasePoiseuilleAnal(gx=gx_air_water, mu_l=1.22*1.5E-05, mu_h=1000 * 1.0E-06, rho_h=1000, rho_l=1.22, h=h)
u_air_water = np.array([pa_air_water.get_u_profile(y_[i]) for i in range(len(y_))])


###################################################################################################################

fig_name = f'two_phase_Poiseuille_anal_vs_fd_rho{rho_h/rho_l}_v{kin_visc_h/kin_visc_l}.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()

plt.plot(u_a, y_,
         color="black", marker="", markevery=5, markersize=5, linestyle="-", linewidth=2,
         label=r'$analytical \, solution$')

plt.plot(u_fd, y_,
         color="black", marker="", markevery=1, markersize=15, linestyle=":", linewidth=2,
         label=r'$FD \, solution$')


# ------ format y axis ------ #
yll = y_.min()
yhl = y_.max()
axes.set_ylim([yll, yhl])
# axes.set_yticks(np.linspace(yll, yhl, 5))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')
# ------ format x axis ------ #
# plt.xlim(0, int(xSIZE / 2))
# plt.xlim(int(xSIZE / 2), xSIZE)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'two phase Poiseuille flow')

plt.xlabel(r'$u_x$')
plt.ylabel(r'$y$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure
