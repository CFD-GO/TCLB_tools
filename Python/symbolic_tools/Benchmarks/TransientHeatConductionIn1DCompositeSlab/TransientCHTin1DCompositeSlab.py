"""
See 'On transient heat conduction in a one-dimensional composite slab.' by Y.Sun and S.Wichman, 2004
"""


from scipy.optimize import newton
from scipy import tan, sqrt, sin, cos, exp
from scipy import array
import numpy as np
import matplotlib.pyplot as plt
import os


class Solver:
    def __init__(self, T1, T0, k, alfas, ds):
        """
        Follows the notation from the article, however the indices starts from 0 to account for python's notation.
        The transient temperature profile across three layer composite slab are to be calculated.
        :param T1: Temperature at the surface of the first layer [K]
        :param T0: Temperature at the surface of the third layer [K]
        :param k: [k0, k1, k2] thermal conductivities of the layers [W/m2K]
        :param alfas: [alfa0, alfa1, alfa2] thermal diffusivities of the layers [m2/s]
        :param ds: [d0, d1, d2] thicknesses of the layers [m]
        """

        self.alfa = alfas
        self.d = ds
        self.t0 = 1  # reference time
        self.k = k
        self.delta = self._calc_deltas()

        self.Delta1 = (self.k[0] / self.k[1]) * sqrt(self.delta[1] / self.delta[0])
        self.Delta2 = (self.k[0] / self.k[2]) * sqrt(self.delta[2] / self.delta[0])

        sum_one_over_ks = (1. / k[0] + 1. / k[1] + 1. / k[2])
        self.DeltaPhi1 = (1. / k[0]) / sum_one_over_ks
        self.DeltaPhi2 = (1. / k[1]) / sum_one_over_ks
        self.DeltaPhi3 = (1. / k[2]) / sum_one_over_ks

        self.T1 = T1
        self.T0 = T0

        self.eigenvalues = []

    def _calc_dzeta(self, x):
        """

        :param x:
        :return: dimensionless distance from the first layer
        """
        dzeta_0 = x/self.d[0]
        dzeta_1 = (x-self.d[0])/self.d[1]
        dzeta_2 = (x - (self.d[0] + self.d[1]))/self.d[2]

        return [dzeta_0, dzeta_1, dzeta_2]

    def _calc_tau(self, t):
        """
        :param t:
        :return: dimensionless time
        """
        return t/self.t0

    def _calc_deltas(self):
        """
        :return: dimensionless thermal diffusivity
        """
        deltas = [self.alfa[i]*self.t0/(self.d[i]*self.d[i]) for i in range(3)]
        return deltas



    def _calc_mi(self, _lamba):
        mi1 = sqrt(self.delta[0] / self.delta[1]) * _lamba
        mi2 = sqrt(self.delta[0] / self.delta[2]) * _lamba
        return mi1, mi2

    def _calc_eigenvalues(self, lamba0 = 1.):
        """
        :param lamba0: starting guess
        :return: eigenvalues from eq.8
        """

        # from scipy.optimize import newton
        # def f(x):
        #     return (1.0 / 4.0) * x ** 3 + (3.0 / 4.0) * x ** 2 - (3.0 / 2.0) * x - 2
        #
        # x = 4
        #
        # x0 = newton(f, x, fprime=None, args=(), tol=1.48e-08, maxiter=50, fprime2=None)
        #
        # print('x: ', x)
        # print('x0: ', x0)
        # print("f(x0) = ", ((1.0 / 4.0) * x0 ** 3 + (3.0 / 4.0) * x0 ** 2 - (3.0 / 2.0) * x0 - 2))

        def f(_lamba):
            # mi1 = sqrt(self.delta[0] / self.delta[1]) * _lamba
            # mi2 = sqrt(self.delta[0] / self.delta[2]) * _lamba

            mi1, mi2 = self._calc_mi(_lamba)

            result = - tan(_lamba) - (self.Delta1*tan(mi1) + self.Delta2*tan(mi2))/(1.0 - (self.Delta2/self.Delta1) * tan(mi1)*tan(mi2))
            return result

        lamba = newton(f, lamba0, fprime=None, tol=1.48e-08, maxiter=50, fprime2=None)
        return lamba

    def _calc_eigenfunctions(self, lamba, dzeta):
        X0 = sin(lamba* dzeta[0])
        alfa = cos(lamba)
        beta = sin(lamba)/self.Delta1

        mi1 = sqrt(self.delta[0] / self.delta[1]) * lamba
        mi2 = sqrt(self.delta[0] / self.delta[2]) * lamba

        X1 = alfa*sin(mi1*dzeta[1]) + beta*cos(mi1*dzeta[1])

        alfa_dash = cos(lamba)*cos(mi1) - sin(lamba)*sin(mi1)/self.Delta1
        beta_dash = cos(lamba)*sin(mi1)*self.Delta1/self.Delta2 + sin(lamba)*cos(mi1)/self.Delta2

        X2 = alfa_dash*sin(mi2*dzeta[2]) + beta_dash*cos(mi2*dzeta[2])

        return [X0, X1, X2]

    def calc_steady_state(self, x):
        """
        :param x: position in physical units [m]
        :return:
        """

        dzeta = self._calc_dzeta(x)
        mask = [0 <= dzeta[i] < 1 for i in range(3)]

        if sum(mask) != 1:
            raise Exception(f"x beyond stab x={x}, dzeta{dzeta}")

        psi = array([1 - self.DeltaPhi1 * dzeta[0],
                     1 - self.DeltaPhi1 - self.DeltaPhi2*dzeta[1],
                     1 - self.DeltaPhi1 - self.DeltaPhi2 - self.DeltaPhi3*dzeta[2]
                     ])  # list indices must be integers or slices, not list

        result = psi[mask][0]  # take the first and only element from the array
        return result

    def remove_duplicates(self, values):
        output = []
        seen = set()
        for value in values:
            # If value has not been encountered yet,
            # add it to both list and set.
            value_rounded = round(value, 8)
            h = hash(value_rounded)
            if h not in seen:
                output.append(value_rounded)
                seen.add(h)
        return output

    def calc_transient_state(self, x, tau):
        dzeta = self._calc_dzeta(x)
        mask = [0 <= dzeta[i] < 1 for i in range(3)]

        if sum(mask) != 1:
            raise Exception(f"x beyond stab x={x}, dzeta{dzeta}")

        N = 10

        initial_guess = 1.4432

        if not self.eigenvalues:
            while len(self.eigenvalues) < N:
                initial_guess += 0.1
                eigenvalue = self._calc_eigenvalues(lamba0=initial_guess)
                self.eigenvalues.append(eigenvalue)
                self.eigenvalues = self.remove_duplicates(self.eigenvalues)

        phi = 0
        for eigenvalue in self.eigenvalues:
            alfa = cos(eigenvalue)
            beta = sin(eigenvalue) / self.Delta1
            X0 = sin(eigenvalue) * dzeta[0]

            mi1, mi2 = self._calc_mi(eigenvalue)
            X1 = alfa * sin(mi1 * dzeta[1]) + beta * cos(mi1 * dzeta[1])

            alfa_dash = cos(eigenvalue) * cos(mi1) - sin(eigenvalue) * sin(mi1) / self.Delta1
            beta_dash = cos(eigenvalue) * sin(mi1) * self.Delta1 / self.Delta2 + sin(eigenvalue) * cos(mi1) / self.Delta2
            X2 = alfa_dash * sin(mi2 * dzeta[2]) + beta_dash * cos(mi2 * dzeta[2])

            c2 = cos(eigenvalue) * cos(eigenvalue)
            s2 = sin(eigenvalue) * sin(eigenvalue)
            s2 /= (self.DeltaPhi1*self.DeltaPhi1)
            M = self.k[1]*self.k[2]/2. \
                + self.k[0]*self.k[2]*(c2 + s2)/2. \
                + self.k[0]*self.k[1]*(alfa_dash*alfa_dash+beta_dash*beta_dash)/2

            A = self.k[1]*self.k[2]/(eigenvalue*M)

            phi_i = array([X0,
                           self.Delta1*X1,
                           self.Delta2*X2])
            eigenfunction = phi_i[mask][0]  # pick the right one

            eigenfunction *= A*exp(-eigenvalue*eigenvalue*self.delta[0]*tau)
            phi += eigenfunction

        return phi

    def calc_dimensionless_temp(self, x, tau):

        psi = self.calc_steady_state(x)
        phi = self.calc_transient_state(x, tau)

        dimensionless_temp = psi + phi
        return dimensionless_temp

    def calc_dimensional_temp(self, x, tau):
        """
        :param x: physical position [m]
        :param tau: nondimensional time argh
        :return: Temperature [K] in position x [m] at
        """
        dimensionless_temp = self.calc_dimensionless_temp(x, tau)
        temp = dimensionless_temp * (self.T1 - self.T0) + self.T0
        return temp


ds = [1., 1., 1.]

solver = Solver(1, 0, k=[1., 0.1, 1.], alfas=[1., 3., 1.], ds=ds)

x = np.linspace(0., 3, num=30, endpoint=False)
x_dimensionless = x
ss = [solver.calc_steady_state(x_i) for x_i in x]

tau = 0.5
ts = [solver.calc_transient_state(x_i, tau=tau) for x_i in x]

if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/TransientHeatConductionInCompositeSlab.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()
plt.plot(x, ss,
         color="black", marker="", markevery=1, markersize=15, linestyle="--", linewidth=2,
         label='analytical solution - steady state')

plt.plot(x, ts,
         color="black", marker="", markevery=1, markersize=15, linestyle=":", linewidth=2,
         label=f'analytical solution - transient contribution at tau={tau}')

# ------ format y axis ------ #
# yll = y.min()
# yhl = y.max()
# axes.set_ylim([yll, yhl])
# axes.set_yticks(np.linspace(yll, yhl, 5))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')


# ------ format x axis ------ #
# plt.xlim(x1-0.5, x2+0.5)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'Sample plot\n '
          # r'$x_{1}$=' + f'{x1}' + '\t' + r'$x_{2}$=' + f'{x2}'
          # f'; \t'
          # r'$x_{step}$' + f'={step:.4f}'
          )
plt.xlabel(r'$x_{label}$')
plt.ylabel(r'$y_{label}$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()


########################################################3
# def f(x, y):
#     return (1.0 / 4.0) * x ** 3 + (3.0 / 4.0) * x ** 2 - (3.0 / 2.0) * x - 2 + y
#
# x = 4
#
# x0 = newton(f, x, fprime=None, args=(1e-5,), tol=1.48e-08, maxiter=50, fprime2=None)
#
# print('x: ', x)
# print('x0: ', x0)
# print("f(x0) = ", ((1.0 / 4.0) * x0 ** 3 + (3.0 / 4.0) * x0 ** 2 - (3.0 / 2.0) * x0 - 2))