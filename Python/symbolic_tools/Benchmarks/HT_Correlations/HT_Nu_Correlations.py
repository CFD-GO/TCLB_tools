import numpy as np

"""
Nusselt number represents the enhancement of heat transfer through a fluid
    as a result of convection relative to conduction across the same fluid layer.
    Nu = h * D / k  = {convective heat transfer }/{conductive heat transfer }
    where h is the average heat transfer coefficient [W/(m2*K)]
    from Newton's law of cooling, h = q/(Area*(T_surf - T_inf))

Prandlt number describes relative thickness of the momentum to thermal boundary layer.
    Pr = v/alfa = mu*cp/k = {molecular diffusivity of momentum BL} / {molecular diffusivity of heat BL}

"""

"""
Source:
'Heat and Mass Transfer' by Yunus A. Cengel, 3rd edition, 2007, McGraw-Hill
Chapter 7 External Forced Convection, p413-414
"""

def get_Nu_cylinder_by_Churchill_Bernstein(Re, Pr):
    """
    For cross flow above a cylinder.
    Valid for Re*Pr > 0.2, smooth surface. best accuracy for Pr ~ 1
    Fluid properties are evaluated at the film temperature T_f =  (T_{surf} + T_{inf})/2
    :param Re:
    :param Pr:
    D is the specific dimension - diameter
    k is heat conductivity [W/(m*K)]
    """
    RePr = Re * Pr
    if np.any(np.less(RePr, 0.2)):
        raise Exception(f"Correlation is valid for Re*Pr > 0.2, smooth surface. RePr={RePr} ")

    Nu_cylinder = (0.62 * pow(Re, 1 / 2) * pow(Pr, 1 / 3))
    Nu_cylinder /= pow((1 + pow(0.4 / Pr, 2 / 3)), 1 / 4)
    Nu_cylinder *= pow(1 + pow(Re / 282000, 5 / 8), 4 / 5)
    Nu_cylinder += 0.3

    return Nu_cylinder

def get_Nu_cylinder_by_Zukauskas_Jacob(Re, Pr):
    if np.any(np.less(Re, 0.4)):
        raise Exception(f"Correlation is not valid for Re < 0.4, smooth surface. Re={Re} ")

    # if np.any(np.greater(Re, 4)):
    #     raise Exception(f"Correlation is valid for Re > 4, smooth surface. Re={Re} ")

    C = 0.989
    m = 0.330
    n = 1/3
    if Re > 4:
        C = 0.911
        m = 0.385
        n = 1 / 3

    if np.any(np.greater(Re, 40)):
        raise Exception(f"Correlation is not valid for Re > 40, smooth surface. Re={Re} ")

    Nu_cylinder = C * pow(Re, m) * pow(Pr, n)
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
    Nu = 0.4 * pow(Re, 2) + 0.06 * pow(Re, 2 / 3)
    Nu *= pow(Pr, 0.4) * pow(mu_inf / mu_surf, 1 / 4)
    Nu += 2
    return Nu
