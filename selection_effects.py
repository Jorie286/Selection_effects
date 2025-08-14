# import necessary packages
import numpy as np
import scipy.constants as const
from scipy.optimize import fsolve
import pygedm
import healpy as hp

# get functions and constants from other files
import gal_cart
import survey
from skytempy import skytemp

# astrophysical constants
R0_Kpc=8.3 # Kpc, distance Sun-Galaxy Center

def DM_fnct(x, y, z, freq):
    """
    Get the dispersion measure of the pulsar.

    Inputs:
    -------
    x, y, z; cartesian coordinates for the position of the pulsar (Kpc)
    freq, observing frequency (MHz)

    Returns:
    --------
    dm, the dispersion measure of the pulsar in the given direction. (pc / cm^3)
    tau_sc, scattering timescale at 1 MHz, units: seconds (from pygedm dist_to_dm docs)
    """

    l, b, d = gal_cart.cart2gal(-x, y, z, degree=True) # Units: l (degrees), b (degrees), d (kpc)

    # d needs to be in pc for the pygedm function
    DM, tau_sc = pygedm.dist_to_dm(l, b, dist = d*1e3, nu=(freq/1e3)) # Units: DM (pc / cm^3), tau_sc (GHz)
    # the values of DM and tau_sc have units attached to them when they are returned by dist_to_dm, remove these
    DM, tau_sc = str(DM).split()[0], str(tau_sc).split()[0]
    return float(DM), tau_sc

def tau_scatt_fnct(DM, freq):
    """
    Get the ISM scattering time for the given pulsar

    Input:
    ------
    DM, dispersion measure of the pulsar (pc/cm^3)
    freq, survey frequency (MHz)

    Returns:
    --------
    tau_scatt, ISM scattering time (seconds)
    """

    tau_scatt = -6.46 + 0.154 * np.log10(DM) + 1.07 * (np.log10(DM) ** 2) - 3.86 * np.log10(freq / 1e3)

    tau_scatt = 10**tau_scatt # ms
    tau_scatt = tau_scatt * 1e-3 # convert to seconds
    return tau_scatt

def DM0_fnct(freq, d_f, n_chan, tau_samp):
    """
    Get the diagonal dispersion measure of the survey.

    Input:
    ------
    freq, survey frequency (MHz)
    d_f, receiver bandwidth (MHz)
    n_chan, number of receiver channels
    tau_samp, sampling time (seconds)

    Returns:
    --------
    DM_0, diagonal dispersion measure (pc cm^-3)
    """
    wavelength=const.c/(freq * 1e6) # Units: m

    DM_0 = 1000.0 * tau_samp *((3e2/wavelength)**3)/(8.3e6*(d_f/n_chan))

    return DM_0

def T_sky_fnct(x, y, z, freq):
    """
    Calculate sky temperature at the pulsar's location.

    Input:
    ------
    x, y, x; cartesian coordinates of the pulsar (Kpc)
    freq, survey frequency (MHz)

    Returns:
    --------
    T_sky, sky temperature in the direction of the puslar (Kelvin)
    """
    l, b, d = gal_cart.cart2gal(x+R0_Kpc, y, z, degree=True) # galactic coordinates of the puslar, Units: l (rad), b (rad), d (kpc)

    # Sky temperature at 408 MHz
    s = skytemp.SkyTemp(l, b, r"./skytempy/haslam408_ds_Remazeilles2014.fits") # get the skytemp information from the fits file
    T_sky = s.get_temp(freq) # get the temperature from the output, freq units in MHz
    return T_sky

def flux(L, x, y, z):
    """
    Compute the flux of the pulsar.

    Input:
        L, luminosity of the modeled pulsar, an input given by user (mJy Kpc^2)
        x, y, z; cartesian coordinates (Kpc)

    Returns:
        F, flux of the pulsar (mJy)
        D, distance to the pulsar (Kpc)
        Area, distance to the pulsar squared (Kpc^2)
    """

    D = np.sqrt(((x-R0_Kpc)**2)+ (y**2) + (z**2)) # distance to pulsar (Kpc)
    Area = D**2 # Kpc^2
    F = L/(4*np.pi*Area) # flux of the pulsar in (mJy)
    return F, D, Area

def S_min(M, P_orb, e, a, P, P_dot, B, x, y, z, vx, vy, vz, L, T_rec, d_f, n_chan, freq, tau_samp, G, t_int, npol = 2, SNmin = 10, beta=1):
    """
    Compute the minimum flux that a pulsar can have and still be detectable.

    Input:
        M, mass (M_sun)
        P_orb, orbital period (days)
        e, eccentricity
        a, separation between the two objects (AU)
        P, rotational period (seconds)
        P_dot, change in rotational period (seconds)
        B, surface magnetic field (Tesla)
        x, y, z; cartesian coordinates (Kpc)
        vx, vy, vz; change in position of each cartesian (km/s)
        L, luminosity (mJy Kpc^2)
        T_rec, receiver temperature (Kelvin)
        d_f, receiver bandwidth (MHz)
        n_chan, number of channels in the survey
        freq, observing frequency (MHz)
        tau_samp, sampling time (seconds)
        G, gain of the telescope (Kelvin/ Jy)
        t_int, integration time (seconds)
        npol, number of polarizations in the detector (automatically set to 2)
        SNmin, minimum detection threshold (automatically set to 10)
        beta, parameter to account for the errors that increase the noise in the signal (automatically set to 1)

    Returns:
        S_min_05, S_min_27, S_min_fwhm; the lower limit of flux a simulated pulsar can have to be detected at a given S/N ratio (mJy) for three different functions of W_i
        F, flux (mJy)
        D, distance (Kpc)
        SNR_05, SNR_27, SNR_fwhm; the signal to noise ratio of the pulsar data for three different functions of W_i
        T_sky, sky temperature (Kelvin)
        DM, dispersion measure (pc cm^-3)
        DM_0, diagonal dispersion measure (pc cm^-3)
        We_05, We_27, We_fwhm; effective pulse width for different intrinsic pulse widths (seconds)
    """

    DM, tau_sc = DM_fnct(x, y, z, freq) # dispersion measure in the direction of the pulsar, Units: pc/cm^3, seconds

    tau_scatt = tau_scatt_fnct(DM, freq) # ISM scattering time, Units: seconds

    DM_0 = DM0_fnct(freq, d_f, n_chan, tau_samp) # diagonal dispersion measure of the survey

    T_sky = T_sky_fnct(-x, -y, z, freq) # get the sky temperature in the direction of the pulsar (Kelvin)

    F, D, Area = flux(L, x, y, z) # get the flux of the pulsar, Units: F (mJy), D (Kpc), Area (Kpc^2)

    # the three equations for W_i come from:
        # fixed duty cycle in DNS paper, list index 0
        # inside C code comment; from https://iopscience.iop.org/article/10.3847/1538-4357/ab75e2/pdf, list index 1
        # from C code, list index 2
    FWHM = 0.04
    Wi_list = [P*0.05, (P**0.27)*0.05, FWHM*P]
    S_min_list = []
    SNR_list = []
    We_list = []

    for W_i in Wi_list:
        # compute the effective pulse width
        We = np.sqrt(W_i**2 + tau_samp**2 + (tau_samp*(DM/DM_0))**2 + tau_scatt**2) # Units: seconds
        We_list.append(We)

        if We>=P:
            S_min = 9.99e9
        else:
            # compute the minimum flux (S_min)
            S_min = beta*((SNmin*(T_rec+T_sky))/((G*1e3)*np.sqrt(npol*t_int*(d_f/1e6))))*np.sqrt(We/(P-We)) # units: mJy

        S_min_list.append(S_min) # add each calculated S_min to the list

        # get the S/N ratio of the pulsar data
        npol = 2
        if We >= P:
            SNR = 0
        else:
            SNR = F / np.sqrt(np.pi / 2) / np.sqrt(We / (P - We)) / (T_rec + T_sky) * (G * np.sqrt(npol * d_f * t_int))
        SNR_list.append(SNR)

    # get each S_min from the list
    S_min_05, S_min_27, S_min_fwhm = S_min_list

    # get each SNR from the list
    SNR_05, SNR_27, SNR_fwhm = SNR_list

    We_05, We_27, We_fwhm = We_list

    return S_min_05, S_min_27, S_min_fwhm, F, Area, SNR_05, SNR_27, SNR_fwhm, T_sky, DM, DM_0, We_05, We_27, We_fwhm

def f_beaming(P):
    """
    Computes the pulsar beaming fraction. It should return a result between 0 and 1.

    Input:
        P, period of the pulsar (seconds)
    Returns:
        f_b, pulsar beaming fraction
    """

    # beaming fraction model, P must be entered in seconds
    f_b = 0.09*((np.log10(P.clip(1e-10)/10))**2)+0.03

    return f_b
