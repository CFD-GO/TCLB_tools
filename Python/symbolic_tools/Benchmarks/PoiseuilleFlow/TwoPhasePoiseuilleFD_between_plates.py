
import numpy as np


def get_tanh_profile(y, h, phi_h, phi_l, W=1):
    """
    :param y: bottom y=-h, top y=h
    :param h: # distance from the center to the channel walls
    :param phi_h: quantity of interest - higher part of the channel
    :param phi_l: quantity of interest - lower part of the channel
    :param W: inteface thickness
    :return:
    """
    y_shifted = y + h  # shift y: bottom y=0, top y=2h
    L = 2*h
    phi = (phi_h+phi_l)/2 + (phi_h-phi_l)/2 * np.tanh((2 * y_shifted - L) / W)
    return phi


class TwoPhasePoiseuilleFD_between_plates:
    """
    CFD lecture notes J. Rokicki, p23
    """
    def __init__(self, gx, mu_l, mu_h, rho_l, rho_h, r):
        self.mu_h = mu_h  # dynamic viscosity of upper fluid
        self.mu_l = mu_l  # dynamic viscosity of lower fluid
        self.rho_h = rho_h  # density of upper fluid
        self.rho_l = rho_l  # density of lower fluid

        self.r = r  # distance from the center to the channel walls
        self.gx = gx  # body force

    def get_u_profile(self, y, W=1):
        """
        :param y: bottom y=-h, top y=h
        :param W: interface width
        :return:
        """

        N = len(y)
        step = 2 * self.r / N
        A = np.identity(N)
        b = np.zeros(shape=(N, 1))

        mu = get_tanh_profile(y, self.r, self.mu_h, self.mu_l, W)
        rho = get_tanh_profile(y, self.r, self.rho_h, self.rho_l, W)

        for i in range(1, N-1):
            b[i] = -self.gx * rho[i] * step * step  # RHS

            A[i, i-1] = mu[i] - (mu[i + 1] - mu[i - 1]) / 4  # left of the diagonal
            A[i, i] = -2 * mu[i]  # the diagonal
            A[i, i+1] = mu[i] + (mu[i+1] - mu[i-1])/4  # right of the diagonal

        b[0] = 0  # first boundary condition
        b[N-1] = 0  # second boundary condition

        u = np.linalg.solve(A, b)
        return u
