# equivalent file of gal_cart.h and gal_cart.c

# import necessary packages and files
import galpy
import galpy.util.coords

def cart2gal(x, y, z):
    """
    Convert the given cartisian coordinates (x, y, z) to galactic coordinates (l, b, d).

    Inputs:
        x, cartisian x coordinate
        y, cartisian y coordinate
        z, cartisian z coordinate

    Returns:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate
        d, galactic d (distance) coordinate
    """

    lbd = galpy.util.coords.XYZ_to_lbd(x, y, z)
    l = lbd[0]
    b = lbd[1]
    d = lbd[2]

    return l, b, d

def gal2cart(l, b, d, x, y, z):
    """
    Convert galactic coordinates (l, b, d) to cartisian coordinates (x, y, z).

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate
        d, galactic d (distance) coordinate

    Returns:
        x, cartisian x coordinate
        y, cartisian y coordinate
        z, cartisian z cooordinate
    """

    xyz = galpy.util.coords.lbd_to_XYZ(l, b, d)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    return x, y, z

def gal_dot(l, b, d, vx, vy, vz):
    """
    Given the galactic coordinates (l, b, d) and the velocities in cartisian coordinates (vx, vy, vz), calculate the velociy in the l and b galactic coordinate (ldot, bdot).

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate
        d, galactic d (distance) coordinate (kpc)
        vx, velocity in cartisian x coordinate (km/s)
        vy, velocity in cartisian y coordinate (km/s)
        vz, velocity in cartisian z coordinate (km/s)

    Returns:
        vr, line of sight velocity (km/s)
        ldot_cosb, change in (velocity of) galactic l coordinate (mas/yr)
        bdot, change in (velocity of) galactic b coordinate (mas/yr)
    """

    r = galpy.util.coords.vxvyvz_to_vrpmllpmbb(vx, vy, vz, l, b, d)
    vr = r[0]
    ldot_cosb = r[1]
    bdot = r[2]

    return r, ldot_cosb, bdot


def gal2eq(l, b):
    """
    Convert galatic coordinates to equatorial coordinates.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        alpha, right ascension angle 0, 2*pi, radians
        delta, delination angle -pi/2, pi/2, radians
    """

    r = galpy.util.coords.lb_to_radec(l, b)
    alpha = r[0]
    delta = r[1]

    return alpha, delta
