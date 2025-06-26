# equivalent file of gal2eq.h and gal2eq.c

# import necessary pacages and files
import numpy as np
import math
import mymath

def gal2eq(l, b):
    """
    Convert galatic coordinates to equatorial coordinates.

    2000.0 ???

    Calculate l, b, distance of pulsar to the Sun. R0_Kpc is the distance Sun-Center of Galaxy, in Kpc, as defined in uconst.py

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        alpha, right ascension angle 0, 2*pi, radians
        delta, delination angle -pi/2, pi/2, radians
    """

    # do we need these?????
    # l_zero = 33 * np.pi / 180
    # alpha_zero = 282.25 * np.pi / 180
    # delta_NGP = (27.4) * np.pi / 180

    # 2000.0 values ???
    l_zero = 32.93 * np.pi / 180
    alpha_zero = 282.86 * np.pi / 180
    delta_NGP = (27 + 7.8 / 60) * np.pi /180

    sin_delta = np.sin(b) * np.sin(delta_NGP) + np.cos(b) * np.cos(delta_NGP) * np.sin(l - l_zero)

    delta = np.arcsin(sin_delta)

    caa = np.cos(l - l_zero) * np.cos(b) / np.cos(delta)
    saa = (-np.sin(b) * np.cos(delta_NGP) + np.cos(b) * np.sin(delta_NGP) * np.sin(l - l_zero)) / np.cos(delta)
    aa = mymath.arg2PI(caa, saa) # ????

    alpha = aa + alpha_zero

    # because of the previous addition, this needs to be reduced again to [0, 2pi)
    alpha = alpha - 2 * np.pi * math.floor(alpha/ (2 * np.pi))

    return alpha, delta
