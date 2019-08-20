"""
See 1.60, 1.61 p34
PhD Thesis `Lattice Boltzmann Simulation of Multiphase Flows;
Application to Wave Breaking and Sea Spray Generation`
by Amir Banari 2014

and
https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_flow_from_the_Navier%E2%80%93Stokes_equations
"""


def calc_gx_between_plates(uc, mu_l, mu_h, rho_l, rho_h, h):
    """
    :param uc: velocity at the interface (center)
    :param mu_l: dynamic viscosity of lower fluid
    :param mu_h: dynamic viscosity of upper fluid
    :param rho_l: density of lower fluid
    :param rho_h: density of upper fluid
    :param h: distance from the center to the channel walls
    :return: gravitation
    """
    gx = 2 * uc * (mu_l + mu_h) / ((rho_l + rho_h) * h * h)
    return gx


def calc_gx_in_pipe(uz_max, rho, nu, D):
    """
    :param uz_max: max velocity at the center of the pipe
    :param nu: kinematic viscosity of the fluid
    :param rho: density of the fluid
    :param D: effective diameter of the pipe
    :return: gravitation
    """
    gx = 16 * uz_max * nu / (rho * D * D)
    return gx


class OnePhasePoiseuilleAnalInPipe:
    def __init__(self, gx, nu, D, rho=1):
        """
        :param gx: gravitation
        :param nu: kinematic viscosity
        :param D: effective diameter of the pipe
        """
        self.gx = gx
        self.nu = nu
        self.R = D/2.
        self.rho = rho

    def get_u_profile(self, r):
        """
        :param r: distance from the center of the pipe
        :return: uz
        """
        uz = self.gx*(self.R*self.R - r*r)/(4*self.nu)
        return uz


class OnePhasePoiseuilleAnalBetweenPlates:
    def __init__(self, gx, nu, H):
        """
        :param gx: gravitation
        :param nu: kinematic viscosity
        :param H: effective height of the channel
        """
        self.gx = gx
        self.nu = nu
        self.H = H

    def get_u_profile(self, y):
        """
        :param y: distance from the bottom wall
        :return: ux
        """
        ux = 0.5 * self.gx * y * (self.H - y) / self.nu
        return ux


class TwoPhasePoiseuilleAnal:
    def __init__(self, gx, mu_l, mu_h, rho_l, rho_h, h):
        self.mu_l = mu_l  # dynamic viscosity of lower fluid
        self.mu_h = mu_h  # dynamic viscosity of upper fluid
        self.rho_l = rho_l  # density of lower fluid
        self.rho_h = rho_h  # density of upper fluid

        self.h = h  # distance from the center to the channel walls
        self.gx = gx  # body force

    def get_u_profile(self, y):
        """
        :param y: distance from the interface (center)
        :return: ux
        """
        if y > 0:
            result = -self.rho_h * (y / self.h) * (y / self.h)
            result -= (y / self.h) * (self.mu_h * self.rho_l - self.mu_l * self.rho_h) / (self.mu_l + self.mu_h)
            result += (self.rho_h + self.rho_l) * self.mu_h / (self.mu_l + self.mu_h)

            result *= self.gx * self.h * self.h / (2 * self.mu_h)

            return result
        else:
            result = -self.rho_l * (y / self.h) * (y / self.h)
            result -= (y / self.h) * (self.mu_h * self.rho_l - self.mu_l * self.rho_h) / (self.mu_l + self.mu_h)
            result += + (self.rho_h + self.rho_l) * self.mu_l / (self.mu_l + self.mu_h)

            result *= self.gx * self.h * self.h / (2 * self.mu_l)
            return result

    def get_uc_from_gx(self, gx=None):

        if gx is None:
            gx = self.gx

        uc = gx * self.h * self.h / (self.mu_l + self.mu_h)
        return uc

