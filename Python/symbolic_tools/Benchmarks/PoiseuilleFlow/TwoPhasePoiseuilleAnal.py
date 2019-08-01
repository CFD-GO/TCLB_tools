"""
See 1.60, 1.61 p34
PhD Thesis `Lattice Boltzmann Simulation of Multiphase Flows;
Application to Wave Breaking and Sea Spray Generation`
by Amir Banari 2014
"""


class OnePhasePoiseuilleAnal:
    def __init__(self, gx, nu, D):
        """
        :param gx: gravitation
        :param nu: kinematic viscosity
        :param D: effective diameter of the channel
        """
        self.gx = gx
        self.nu = nu
        self.D = D

    def get_u_profile(self, y):
        """
        :param y: distance from the bottom wall
        :return: ux
        """
        ux = 0.5*self.gx*y*(self.D - y)/self.nu
        return ux


def calc_gx(uc, mu_l, mu_h, rho_l, rho_h, r):
    """
    :param uc: velocity at the interface (center)
    :param mu_l: dynamic viscosity of lower fluid
    :param mu_h: dynamic viscosity of upper fluid
    :param rho_l: density of lower fluid
    :param rho_h: density of upper fluid
    :param r: distance from the center to the channel walls
    :return: gravitation
    """
    gx = 2 * uc * (mu_l + mu_h) / ((rho_l + rho_h) * r * r)
    return gx


class TwoPhasePoiseuilleAnal:
    def __init__(self, gx, mu_l, mu_h, rho_l, rho_h, r):
        self.mu_l = mu_l  # dynamic viscosity of lower fluid
        self.mu_h = mu_h  # dynamic viscosity of upper fluid
        self.rho_l = rho_l  # density of lower fluid
        self.rho_h = rho_h  # density of upper fluid

        self.r = r  # distance from the center to the channel walls
        self.gx = gx  # body force

    def get_u_profile(self, y):
        """
        :param y: distance from the interface (center)
        :return: ux
        """
        if y > 0:
            result = -self.rho_h * (y / self.r) * (y / self.r)
            result -= (y / self.r) * (self.mu_h * self.rho_l - self.mu_l * self.rho_h) / (self.mu_l + self.mu_h)
            result += (self.rho_h + self.rho_l) * self.mu_h / (self.mu_l + self.mu_h)

            result *= self.gx * self.r * self.r / (2 * self.mu_h)

            return result
        else:
            result = -self.rho_l * (y / self.r) * (y / self.r)
            result -= (y / self.r) * (self.mu_h * self.rho_l - self.mu_l * self.rho_h) / (self.mu_l + self.mu_h)
            result += + (self.rho_h + self.rho_l) * self.mu_l / (self.mu_l + self.mu_h)

            result *= self.gx * self.r * self.r / (2 * self.mu_l)
            return result

    def get_uc_from_gx(self, gx=None):

        if gx is None:
            gx = self.gx

        uc = gx * self.r * self.r / (self.mu_l + self.mu_h)
        return uc

