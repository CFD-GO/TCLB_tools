#  See chapter 8.6.2, eq 8.67, p324 from
# 'The Lattice Boltzmann Method: Principles and Practice'
#  by T.Krüger, H.Kusumaatmaja, A.Kuzmin, O.Shardt, G.Silva, E.M.Viggen
# which origins from section 2.1 p334 from
# 'Hydrodynamics, Mass and Meat Transfer in Chemical Engineering'
# A.D. Polyanin, A.M. Kutepov, A.V. Vyazmin and D.A. Kazenin


# 2D, unsteady state case:
# Unsteady heat conduction inside a cylinder.
# The cylinder is initialized with initial temperature = 0.
# Dirichlet BC is prescribed on the wall of the cylinder and the cylinder warms up.



import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
# import scipy.special.j0 as J0
# import scipy.special.j1 as J1

# it is assumed that the initial concentration is 0
cc = 1  # concentration at the boundary
a = 40  # [lu] radius of the cylinder
D = 0.0052  # [lu2/ts] diffusivity
normalized_time = 0.1095  # normalized_time = D * t / (a * a)
normalized_time2 = 0.0547  # second plot

# [2.4048, 5.5201, 8.6537, 11.7915, 14.9309] # first five roots of the zeroth order Bessel function
# as n --> inf the difference between subsequent roots is pi.
n_bessel_roots = 10
bessel_roots = special.jn_zeros(0, n_bessel_roots)


def concentration(r, _normalized_time):
    temp = 0
    for i in range(n_bessel_roots):
        temp += 2. * special.j0(bessel_roots[i] * r / a) \
                * np.exp(-bessel_roots[i] * bessel_roots[i] * _normalized_time) \
                / (bessel_roots[i] * special.j1(bessel_roots[i]))

    # for i in range(n_bessel_roots):
    #     temp += 2. * special.spherical_jn(0, bessel_roots[i] * r / a) \
    #             * np.exp(-bessel_roots[i] * bessel_roots[i] * D * t / (a * a))\
    #             / (bessel_roots[i] * special.spherical_jn(1, bessel_roots[i]))

    result = cc*(1-temp)
    return result


x = np.linspace(0, a, 100)
y = [concentration(i, normalized_time) for i in x]
y2 = [concentration(i, normalized_time2) for i in x]

# make plot
plt.rcParams.update({'font.size': 18})
plt.figure(figsize=(14, 9))

plt.plot(x / a, y, color="black", linestyle=":", linewidth=3, label='Analytical $\hat{t}$=' + f'{normalized_time}')
plt.plot(x / a, y2, color="black", marker=">", markersize=9, linestyle="",
         label='Analytical $\hat{t}$=' + f'{normalized_time2}')

# plt.plot(u_fd, y_fd, color="black", linestyle="-",  linewidth=2, label='FD - diffuse interface')
# plt.plot(u_diff, (x_diff - len(x_diff)/2 + 0.5), color="black", marker="o", markersize=9, linestyle="", label='current model - diffuse interface')

axes = plt.gca()
axes.set_yticks(np.arange(0, 1.0001, 0.2))
axes.set_ylim([0, 1])
plt.xlim(0, 1)
axes.set_xticks(np.arange(0, 1.0001, 0.2))

#     plt.ylim(y1.min(), y1.max())
#     plt.ylim(1.25*min(y1.min(), y2.min()), 1.25*max(y1.max(), y2.max()))

plt.ylabel(r'$C/C_c$')
plt.xlabel(r'$r/a$')

plt.title('Diffusion from cylinder without flow')
plt.grid(True)
plt.legend()

# plt.text(0.0, 5E-6, r'$\mu^* = %s$' % str(mu_ratio))
# plt.text(0.0, 5E-6, r'$\rho^* = %s$' % str(rho_l/rho_g))
fig = plt.gcf()  # get current figure
fig.savefig(f'Diffusion_from_cylinder_without_flow_t1{normalized_time}_t2{normalized_time2}.pdf', bbox_inches='tight')
plt.show()
