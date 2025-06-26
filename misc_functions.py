# equivalent of misc_functions.h and misc_functions.c

import numpy as np
import math
import fortran_functions

def print_usage():
    printf("\nUsage:\n\n")
    printf("psrevolve [<parameter1> <parameter2> ... <parameter20>]\n\n")
    printf("The 20 parameters required are as follows:\n\n")
    printf("   1: Number of sets of pulsars (1000 pulsars each) to be processed.\n")
    printf("   2: logP_mean,  dec. log. of mean value of pulsar period P  [ log10 s ]\n")
    printf("   3: logP_sigma, dec. log. of stdandard deviation of P values [ log10 s ]\n")
    printf("   4: logB_mean,  dec. log. of mean value of pulsar magnetic field B [ log10 G ]\n")
    printf("   5: logB_sigma, dec. log. of standard deviation of B values        [ log10 G ]\n")
    printf("   6: rsl, radial scale length [ Kpc ]\n")
    printf("   7: zsl, z-axis scale length [ Kpc ]\n")
    printf("   8: kick_1d, one dimensional velocity standard deviation [ Km/s ]\n")
    printf("   9: loggamma_mean,  dec. log. of mean val. of pulsar gamma factor [log10 mJy Kpc^2 ]\n")
    printf("  10: loggamma_sigma, dec. log. of std. dev. of gamma factor values [log10 mJy Kpc^2 ]\n")
    printf("  11: sp_idx_mean,  mean value of spectral index of pulsar [ dimensionless ]\n")
    printf("  12: sp_idx_sigma, std.dev. of spectral index values      [ dimensionless ]\n")
    printf("  13: age_min, minimum value of pulsar age [ Myr ]\n")
    printf("  14: age_max, maximum value of pulsar age [ Myr ]\n")
    printf("  15: alpha_mean, mean value of period exponent alpha   [ val.dep.units ]\n")
    printf("  16: alpha_sigma, std. dev. of alpha values            [ val.dep.units ]\n")
    printf("  17: beta_mean, mean value of Pdot exponent beta [ dimensionless ]\n")
    printf("  18: beta_sigma, std. dev. of beta values        [ dimensionless ]\n")
    printf("  19: FWHM_min, min. val of full width at half maximum [ dimensionless ]\n")
    printf("  20: FWHM_max, max. val of full width at half maximum [ dimensionless ]\n")
    printf("\n")

def vcirc(vc_r, vc_z):
    """
    dpdrout is in kpc. solar_masses / (Myr.Myr)

    Mass part goes away as mass is 1 solar mass

    Returns value in km/sec
    """

    vc_dpdr, vc_dpdz = fortran_functions.kggalmod(vc_r, vc_z)
    return (np.sqrt(vc_dpdr * vc_r) * 977.79235)
