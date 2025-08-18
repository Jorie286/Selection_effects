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

def DNS_NSBH_sel_eff(P1, P2, P_orb, e, type1, type2):
    """
    Get the radio detectability cutoff for (Double Neutron Star) DNS and (Neutron Star Black Hole) NSBH systems using the parameters described in Chattopadhyay, D., Stevenson, S., Hurley, J. R., Rossi, L. J., & Flynn, C. 2020, Monthly Notices of the Royal Astronomical Society, 494, 1587, doi: 10.1093/mnras/staa756

    Input:
    ------
    P1, spin period of primary object (seconds)
    P2, spin period of secondary object (seconds)
    P_orb1, orbital period (days)
    e, eccentricity of the binary systems
    type1, the type of the first object
    type2, the type of the second object

    Returns:
    --------
    alt_det1, binary variable indicating if the primary object is radio detectable
    alt_det2, binary varible indictating if the secondary object is radio detectable
    """

    # check to see that neither object in the system is white dwarf; if not, compute radio detectability, else return None for both values
    if type1 != 10 or type1 != 11 or type1 != 12 or type2 != 10 or type2 != 11 or type2 != 12 or type1.lower().find("wd") == -1 or type2.lower().find("wd") == -1:

        # check for a double neutron star system
        if (type1 == 13 and type2 == 13) or (type1.lower().find("ns") != -1 and type2.lower().find("ns") != -1):

            # constants from linear regression fitting with both neutron stars of mass 1.4 M_sun, 1000 s observations and 60 degree orbital inclination
            m_m = -8.90
            c_m= -27.68
            m_c= -3.40
            c_c= 5.72

            m = (m_m * e) + c_m
            c = (m_c * e) + c_c

            alt_det1 = (P_orb >= (m * P1) + c).astype(int)
            alt_det2 = (P_orb >= (m * P2) + c).astype(int)


        else: # the only other system possible is a neutron star black hole system, compute if not a white dwarf neutron star system or a double neutron star system

            # constants from linear regression fitting with black hole mass of 10 M_sun and neutron stars of mass 1.4 M_sun, 1000 s observations and 60 degree orbital inclination
            m_m = -26.42
            c_m= -18.31
            m_c= -2.53
            c_c= 4.51

            m = (m_m * e) + c_m
            c = (m_c * e) + c_c

            alt_det1 = (P_orb >= (m * P1) + c).astype(int)
            alt_det2 = (P_orb >= (m * P2) + c).astype(int)

    else: # if the system does contain a white dwarf, don't get the radio detectability
        alt_det1=None
        alt_det2=None

    return alt_det1, alt_det2

def death_lines(P, P_dot, freq, area, L):
    """
    Determine if the pulsar falls within or outside of the death lines.

    Input:
    ------
    P, spin period (s)
    P_dot, change in the spin period (s/s)
    freq, observing frequency (MHz)
    area, reciever area (m^2) ????
    L, radio luminosity (mJy)

    Returns:
    --------
    P, original spin period if above death lines or None if below death lines
    """

    # determine if the
    death_line = (10**(0.92 * np.log10(P) - 18.65)) >= P_dot

    # DEBATRI: What is the moment of interia (I) here?
    I = 1e38 # kg m^2 (from Lorimer, D., et. al., Handbook of Pulsar Astronomy)
    E_dot = 4 * (np.pi**2) * I * P_dot * (P**(-3)) # Watts
    E_dot = 1e26 * (E_dot/(area * (1e6*freq))) # convert E_dot to Jy
    E_dot = E_dot*1e3 # convert E_dot to mJy to match with L


    # get the radio efficiency of the pulsars
    xi = L/E_dot

    xi_max = 0.01

    # check the death line and radio efficieny of the pulsar, if it has stopped emitting, set P_1 to None
    if death_line == True or xi >= xi_max:
        P = None

    # Return the changed (or unchanged) value of P so we know what pulsar are below the death lines, don't calculate selection effects for these
    return P

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
        P, spin period (seconds)
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

    for W_i in Wi_list:
        # compute the effective pulse width
        We = np.sqrt(W_i**2 + tau_samp**2 + (tau_samp*(DM/DM_0))**2 + tau_scatt**2) # Units: seconds

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

    return S_min_05, S_min_27, S_min_fwhm, F, Area, SNR_05, SNR_27, SNR_fwhm, T_sky

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
