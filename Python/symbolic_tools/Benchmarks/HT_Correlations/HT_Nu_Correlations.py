"""
Nusselt number represents the enhancement of heat transfer through a fluid
    as a result of convection relative to conduction across the same fluid layer.
    Nu = h * D / k  = {convective heat transfer }/{conductive heat transfer }
    where h is the average heat transfer coefficient [W/(m2*K)]
    from Newton's law of cooling, h = q/(Area*(T_surf - T_inf))

Prandlt number describes relative thickness of the momentum to thermal boundary layer.
    Pr = v/alfa = mu*cp/k = {molecular diffusivity of momentum BL} / {molecular diffusivity of heat BL}

"""

def get_Nu_cylinder_by_Churchill_Bernstein(Re, Pr):
    """
    For cross flow above a cylinder.
    Valid for Re*Pr > 0.2, smooth surface.
    Fluid properties are evaluated at the film temperature T_f =  (T_{surf} + T_{inf})/2
    :param Re:
    :param Pr:
    D is the specific dimension - diameter
    k is heat conductivity [W/(m*K)]
    """
    RePr = Re*Pr
    if RePr < 0.2:
        raise Exception(f"Correlation is valid for Re*Pr > 0.2, smooth surface. RePr={RePr} ")

    Nu_cylinder = (0.62 * pow(Re, 2) * pow(Pr, 3))/pow((1+pow(0.4*Pr, 2/3)), 1/4)
    Nu_cylinder *= pow(1 + pow(Re/282000, 5/8), 4/5)
    Nu_cylinder += 0.3

    return Nu_cylinder


def get_Nu_sphere_by_Whitaker(Re, Pr, mu_inf, mu_surf):
    """
    Valid for 3.5 < Re < 80,000 and 0.7 < Pr < 380, smooth surface.
    The fluid properties are evaluated at the free‐stream temperature T∞,
    except for μs which is evaluated at surface temperature.
    :param Re:
    :param Pr:
    :param mu_inf:
    :param mu_surf:
    :return:
    """
    Nu = 0.4*pow(Re, 2) + 0.06*pow(Re, 2/3)
    Nu *= pow(Pr, 0.4)*pow(mu_inf / mu_surf, 1 / 4)
    Nu += 2
    return Nu