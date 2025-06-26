# equivalent file of gal_cart.h and gal_cart.c

# import necessary packages and files
import numpy as np
import math
import uconst

#Calculate l, b, distance of pulsar to the Sun. R0_Kpc is the distance Sun-Center of Galaxy, in Kpc, as defined in uconst.py

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

    # printf("cart2gal called with x=%3.3f, y=%3.3f, z=%3.3f\n", %(x, y, z))

    dx = uconst.R0_Kpc - x
    dxy = np.sqrt(dx * dx + y *y)

    l = 180.0 / np.pi * float(np.arctan2((-1.0)*y, dx))

    b = 180.0 / np.pi * float(np.arctan2(z, dxy))

    d = np.sqrt(dx * dx + y * y + z * z)


    l_before = l
    # reduce l to [0, 360)
    l = l - 360.0 * math.floor(l/360) # check !!!!

    # if l == 360:
    #     printf("l = 360 after reducing it!!! l_before=%3.3f\n", l_before)
    # while l = -0.000008:
    #     print("We get 360!!!")
    # We may not need the following line!
    if (l==360):
        l=0

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

    R0 = 8.5

    # check l and b
    if l < 0 or l >= 360 or math.fabs(b) > 90:
        printf("ERROR: gal2cart():\n")
        printf("ERROR: l=%3.3f, b=%3.3f\n", % (l, b))
        printf("ERROR: l should be inside [0, 360) deg, b should be inside [-90, 90] deg.\n")
        break

    x = uconst.R0_Kpc - d * np.cos(b * np.pi / 180) * np.cos(l * np.pi / 180)
    y = -d * np.cos(b * np.pi / 180) * np.sin(l * np.pi / 180)
    z = d * np.sin(b * np.pi / 180)

    return x, y, z

def gal_dot(l, b, d, vx, vy, vz, ldot, bdot):
    """
    Given the galactic coordinates (l, b, d) and the velocities in cartisian coordinates (vx, vy, vz), calculate the velociy in the l and b galactic coordinate (ldot, bdot).

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate
        d, galactic d (distance) coordinate
        vx, velocity in cartisian x coordinate
        vy, velocity in cartisian y coordinate
        vz, velocity in cartisian z coordinate

    Returns:
        ldot, change in (velocity of) galactic l coordinate
        bdot, change in (velocity of) galactic b coordinate
    """

    # check l and b
    if l < 0 or l >= 360 or fabs(b) > 90: # Check!!!!
        printf("ERROR: gal2cart():\n")
        printf("ERROR: l=%3.3f, b=%3.3f\n", % (l, b))
        printf("ERROR: l should be inside [0, 360) deg, b should be inside [-90, 90] deg.\n")
        break

    factor = (21.0 / 100.0) / d
    ldot = factor * (-np.sin(l * np.pi / 180) * vx + np.cos(l* np.pi / 180) * vy)
    bdot = factor * (-np.sin(b * np.pi / 180) * np.cos(l* np.pi / 180) * vx - np.sin(b * np.pi / 180) * np.sin(l * np.pi / 180) * vy + np.cos(b * np.pi / 180) * vz)

    return ldot, bdot
