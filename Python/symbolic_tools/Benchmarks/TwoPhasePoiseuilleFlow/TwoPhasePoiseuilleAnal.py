"""
See 1.60, 1.61 p34
PhD Thesis `Lattice Boltzmann Simulation of Multiphase Flows;
Application to Wave Breaking and Sea Spray Generation`
by Amir Banari 2014
"""


def calc_gx(uc, mu_l, mu_h, rho_l, rho_h, h):
    """
    :param uc: velocity at the interface (center)
    :param mu_l: dynamic viscosity of lower fluid
    :param mu_h: dynamic viscosity of upper fluid
    :param rho_l: density of lower fluid
    :param rho_h: density of upper fluid
    :param h: distance from the center to the channel walls
    :return:
    """
    gx = 2 * uc * (mu_l + mu_h) / ((rho_l + rho_h) * h * h)
    return gx


class TwoPhasePoiseuilleAnal:
    def __init__(self, gx, mu_l, mu_h, rho_l, rho_h, h):
        self.mu_l = mu_l  # dynamic viscosity of lower fluid
        self.mu_h = mu_h  # dynamic viscosity of upper fluid
        self.rho_l = rho_l  # density of lower fluid
        self.rho_h = rho_h  # density of upper fluid

        self.h = h  # distance from the center to the channel walls
        self.gx = gx  # body force

    def get_u_profile(self, y):
        if y > 0:
            result = -self.rho_h * (y / self.h) * (y / self.h)
            result -= (y / self.h) * (self.mu_h * self.rho_l - self.mu_l * self.rho_h) / (self.mu_l + self.mu_h)
            result += (self.rho_h + self.rho_l) * self.mu_h / (self.mu_l + self.mu_h)

            result *= self.gx * self.h * self.h / (2 * self.mu_h)

            return result
        else:
            result = -self.rho_l*(y / self.h) * (y / self.h)
            result -= (y / self.h) * (self.mu_h * self.rho_l - self.mu_l * self.rho_h) / (self.mu_l + self.mu_h)
            result += + (self.rho_h + self.rho_l) * self.mu_l / (self.mu_l + self.mu_h)

            result *= self.gx * self.h * self.h / (2 * self.mu_l)
            return result

    def get_uc_from_gx(self, gx=None):

        if gx is None:
            gx = self.gx

        uc = gx * self.h * self.h / (self.mu_l + self.mu_h)
        return uc

