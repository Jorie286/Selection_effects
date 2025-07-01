# equivalent file of survey.h and survey.c

import numpy as np
import math
import usconst
import gal2eq

# define an array that includes all of the details of the surveys
# Field name: name    Trec    recBW(MHz)    n_chan    freq(MHz)    t_samp(s)    gain    t_int(s)    is_covered()
surv_array=np.array(["Parkes_70_cm",      60,   32,   256,   436,        0.0003,    0.65, 157,      ic_parkes_70],
                    ["Parkes_Multi_Beam", 24,   288,  96,    1374,       0.00025,   0.65, 2100,     ic_parkes_mb,
                    ["Swin_Interm Lat",   24,   288,  96,    1374,       0.000125,  0.65, 265,      ic_swin_il],
                    ["Swin_Extended",     24,   288,  96,    1374,       0.000125,  0.65, 265,      ic_swin_ext],
                    ["Burgay_et_al",      24,   288,  96,    1374,       0.000125,  0.65, 265,      ic_burgay],
                    ["GMRT",              90,   32,   512,    610,       0.000125,  0.65*15*(45.0/64.0*45.0/64.0), 300, ic_gmrt],
                    ["Parkes_Multi_Beam_ALLSKY", 24,   288,  96,    1374,       0.00025,   0.65, 2100,     ic_all_sky],
                    ["Parkes_Multi_Beam_part", 24,   288,  96,    1374,       0.00025,   0.65, 2100,     ic_part_sky],
                    ["MeerKat",           18, 800, 1024, 1400, 0.000064, 2.8, 2100, ic_all_sky_MeerSKA],
                    ["SKA",               30, 300, 1024, 1400, 0.000064, 8.4, 2100, ic_all_sky_MeerSKA],
                    ["MeerKat_GalPl", 18, 800, 1024, 1400, 0.000064, 2.8, 2100, ic_GalPlane_MeerSKA],
                    ["MeerKat_tint", 18, 800, 1024, 1400, 0.000064, 2.8, 300, ic_all_sky_MeerSKA],
                    ["MeerKat_Galpl_tint", 18, 800, 1024, 1400, 0.000064, 2.8, 300, ic_GalPlane_MeerSKA],
                    ["Lowlat", 23, 340, 870,1240 ,0.000256, 0.65, 4320, ic_lowlat],
                    ["TRUMP_Meer", 18, 400, 1024, 1400, 0.000064, 1.83, 637, ic_TRUMP_Meer],
                    ["MeerKat", 18, 800, 1024, 1400, 0.000064, 2.8, 2100, ic_all_sky_MeerSKA ],
                    ["SKA",     30, 300, 1024, 1400, 0.000064, 8.4, 2100, ic_all_sky_MeerSKA ],
                    ["MeerKat", 18, 800, 1024, 1400, 0.000064, 2.8, 300, ic_GalPlane_MeerSKA ],
                    ["SKA",     30, 300, 1024, 1400, 0.000064, 8.4, 300, ic_GalPlane_MeerSKA ])

# Geometry checking functions  - for each survey
def ic_parkes_70(l, b):
    """
    Parkes 70 cm geometry check.

    Input:
        l, galactic l coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """
    double alpha, delta;

    # This survey is defined in equatorial coordinates.
    alpha, delta = gal_cart.gal2eq(l*np.pi/180, b*np.pi/180)
    alpha=alpha*12/np.pi
    delta=delta*180/np.pi
    # printf("alpha=%.3lf hours, delta=%.3lf deg\n", alpha,delta)

    if ( delta <=  0   and   0 <= alpha and alpha <= 24 ):
        return 1
    else:
        return 0


def ic_parkes_mb(l, b):
    """
    Parkes Multi Beam geometry check.

    Input:
        l, galactic l coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """

    if (math.fabs(b)<=5 and 0 <=l and l <= 50) or (260 <= l and l <360):
        return 1
    else:
        return 0

def ic_swin_il(l, b):
    """
    Swin Interm Lat geometry check.

    Input:
        l, galactic l coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """
   if ( 5 <= math.fabs(b) and math.fabs(b) <= 15 and (  (0 <= l and l <= 50) or (260 <=l and l < 360) ) ):
      return 1
   else:
      return 0

def ic_swin_ext(l, b):
    """
    Swin Extended geometry check.

    Input:
        l, galactic l coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """
    if (  ( 15 <= math.fabs(b) and math.fabs(b) <= 30 and (  (0 <= l and l <= 50) or (260 <=l and l < 360) ) )):
      return 1
   else:
      return 0

def ic_burgay(l, b):
    """
    Burgay et. al. geometry check.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """
    if( math.fabs(b) <= 60   and  220 <= l and l <= 260  ):
       return 1
    else:
       return 0

def ic_gmrt(l, b):
    """
    GMRT geometry check.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """
    # This survey is defined in equatorial coordinates.
    alpha, delta = gal_cart.gal2eq(l*np.pi/180, b*np.pi/180)
    alpha=alpha*12/np.pi
    delta=delta*180/np.pi
    # printf("alpha=%.3lf hours, delta=%.3lf deg\n", alpha,delta)

    if ( delta > -40   and  0 <= alpha and alpha <= 24 ):
       return 1
    else:
       return 0

def ic_all_sky(l, b):
    """
    Parkes Multi Beam ALLSKY geometry check.
    """
    return 1

def ic_part_sky(l, b):
    """
    Parkes Multi Beam part geometry check.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """
    if (math.fabs(b) <= 30 and (0 <= l and l <= 90 or 200 <= l and l <= 360) ):
        return 1
    else:
        return 0

def ic_all_sky_MeerSKA(l, b)
    """
    MeerKat tint geometry check.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """
    # This survey is defined in equatorial coordinates.
    alpha, delta = gal_cart.gal2eq(np.pi/180, b*np.pi/180)
    alpha=alpha*12/np.pi
    delta=delta*180/np.pi

    if ( delta < 30 ):
       return 1
    else:
       return 0

def ic_GalPlane_MeerSKA(l, b):
    """
    MeerKat Galpl tint geometry check.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """

    # This survey is defined in equatorial coordinates.
    alpha, delta = gal_cart.gal2eq(l*np.pi/180, b*np.pi/180)
    alpha=alpha*12/np.pi
    delta=delta*180/np.pi

    if ( delta < 5 and delta > -5 ):
       return 1
    else:
       return 0

def ic_TRUMP_Meer(l, b):
    """
    TRUMP Meer geometry check.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """

    if (math.fabs(b) <= 5.2 and -110 <= l and l <= 10):
       return 1
    else:
       return 0


def ic_lowlat(l, b):
    """
    Lowlat geometry check.

    Input:
        l, galactic l (longitude) coordinate
        b, galactic b coordinate

    Returns:
        1 or 0 (statement is True or False)
    """

    if (math.fabs(b) <= 3.5 and -80 <= l and l <= 30):
       return 1
    else:
       return 0
