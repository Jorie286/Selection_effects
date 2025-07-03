# equivalent file of pulsar.h and pulsar.c

import numpy as np
import math
import misc_functions
import fortran_functions
import gal_cart

def print_pulsar(p):
    """
    Prints the values of all the pulsar member variables

    Input:
        p, ????
    """

    printf("-----------------------\n")
    printf("PULSAR member variables\n")
    printf("-----------------------\n")

    # Always use %g - we never know what the value may be (there may be bugs).
    printf("(x,y,z)=(%g, %g, %g)\n", p.x, p.y, p.z)
    printf("(vx,vy,vz)=(%g, %g, %g)\n", p.vx, p.vy, p.vz)
    printf("age=%g\n", p.age)
    printf("P=%g\n",p.P)
    printf("B=%g\n", p.B)
    printf("gamma_factor=%g\n", p.gamma_factor)
    printf("alpha=%g, beta=%g\n", p.alpha, p.beta)
    printf("spectral index=%g\n", p.sp_idx)
    printf("S400=%g\n", S400_fnct(p))
    printf("DM=%g\n", DM_fnct(p))
    # printf("W_i=%g\n", p.W_i)
    printf("FWHM=%g\n", p.FWHM)
    printf("bmask=%d\n", p.bmask)
    printf("(x0,y0,z0)=(%g, %g, %g)\n", p.x0, p.y0, p.z0)
    printf("(vx0,vy0,vz0)=(%g, %g, %g)\n", p.vx0, p.vy0, p.vz0)
    printf("P0=%g\n", p.P0)
    printf("B0=%g\n", p.B0)
    printf("Calculated variables:\n")
    gal_cart.cart2gal(p.x, p.y, p.z)
    printf("d=%g\n", d)
    printf("L400_fnct(p)=%g\n", L400_fnct(p))
    printf("gal_dot(p)=%g\n", gal_dot(p))

    printf("\n")

def read_pset_cmdline(argc, arg_list):
    """
    Reads the parameter set from command line.

    Input:
        argc, number of parameters
        arg_list, list of parameters

    Returns:
        pset, set of parameters
    """

def print_pset(pset):
    printf("logP_mean=%g\n", pset.logP_mean)
    printf("logP_sigma=%g\n", pset.logP_sigma)

    printf("logB_mean=%g  [log G]\n", pset.logB_mean)
    printf("logB_sigma=%g [log G]\n", pset.logB_sigma)

    printf("rsl=%g\n", pset.rsl)
    printf("zsl=%g\n", pset.zsl)

    printf("kick_1d=%g\n", pset.kick_1d)

    printf("loggamma_mean=%g\n", pset.loggamma_mean)
    printf("loggamma_sigma=%g\n", pset.loggamma_sigma)

    printf("sp_idx_mean=%g\n", pset.sp_idx_mean)
    printf("sp_idx_sigma=%g\n", pset.sp_idx_sigma)

    printf("age_min=%g [Myr]\n", pset.age_min)
    printf("age_max=%g [Myr]\n", pset.age_max)

    printf("alpha_mean=%g\n", pset.alpha_mean)
    printf("alpha_sigma=%g\n", pset.alpha_sigma)

    printf("beta_mean=%g\n", pset.beta_mean)
    printf("beta_sigma=%g\n", pset.beta_sigma)

    printf("FWHM_min=%g\n", pset.FWHM_min)
    printf("FWHM_max=%g\n", pset.FWHM_max)

def print_result(pr, result):
