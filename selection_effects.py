# radio selection effects functions

# import necessary packages
import numpy as np
import scipy.constants as const
import pygedm

# get functions and constants from other files
import gal_cart
import survey
from skytempy import skytemp

# astrophysical constants
R0_Kpc=8.5 # Kpc, distance Sun-Galaxy Center

def DM_fnct(x, y, z):
    """
    Get the dispersion measure of the pulsar.

    Inputs:
    -------
    x, y, z; cartesian coordinates for the position of the pulsar
    freq, observing frequency (MHz)

    Returns:
    --------
    dm, the dispersion measure of the pulsar in the given direction. (pc / cm^3)
    tau_sc, scattering timescale at 1 MHz (from pygedm dist_to_dm docs)
    """

    l, b, d = gal_cart.cart2gal(x, y, z) # Units: l (rad), b (rad), d (kpc)
    d = d * 1e3 # convert d from kpc to pc

    # NOTE: d needs to be in pc for the pygedm function
    DM, tau_sc = pygedm.dist_to_dm((l, b), dist = d, nu=(freq/1e3)) # Units: DM (pc / cm^3), tau_sc (GHz)
    tau_sc = tau_sc*1e3 # convert tau_sc to MHz timescale

    return DM, tau_sc

def tau_scatt_fnct(DM, freq):
    """
    Get the ISM scattering time for the given pulsar

    Input:
    ------
    DM, dispersion measure of the pulsar
    freq, survey frequency (MHz)

    Returns:
    tau_scatt, ISM scattering time (seconds)
    """
    # NOTE: need to check what units need to go in here since unit analysis will not work with log10.
    tau_scatt = -6.46 + 0.154 * np.log10(DM) + 1.07 * np.log10(DM) ** 2 - 3.86 * np.log10(freq / 1e3)
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
    DM_0, diagonal dispersion measure
    """
    wavelength=const.c/(freq * 1e6) # m
    # printf("wavelength=%3.3f\n", % wavelength)
    DM_0 = 1000.0 *tau_samp*((3e2/wavelength)**3)/(8.3e6*(d_f/n_chan)) # NOTE: units ??? s^2 m^-3 ???, should be m^-3 ????
    return DM_0

def T_sky_fnct(x, y, z, freq):
    """
    Calculate sky temperature at the pulsar's location.

    Input:
    ------
    x, y, x; cartesian coordinates of the pulsar
    freq, survey frequency (MHz)

    Returns:
    --------
    T_sky, sky temperature in the direction of the puslar
    """
    l, b, d = gal_cart.cart2gal(x, y, z) # galactic coordinates of the puslar, Units: l (rad), b (rad), d (kpc)
    # Sky temperature at 408 MHz
    s = skytemp.SkyTemp(l, b, '.\\skytempy\\haslam408_ds_Remazeilles2014.fits') # get the skytemp information from the fits file
    T_sky = s.get_temp(freq) # get the temperature from the output, freq units in MHz

    return T_sky

def flux(L, x, y, z):
    """
    Compute the flux of the pulsar.

    Input:
        L, luminosity of the modeled pulsar, an input given by user (Watts ???)
        x, y, z; cartesian coordinates

    Returns:
        F, flux of the pulsar (mJy)
        D, distance to the pulsar (Kpc)
    """

    D = x-R0_Kpc # distance to pulsar (Kpc)
    F = L/(4*np.pi*(D**2)) # flux of the pulsar in (Watts / Kpc^2) ????
    return F, D

def eccentricity(P, a, e, i, T, . . . ):
    """
    Compute the effects eccentricity in a binary orbit has on the detectability of a pulsar. A pulsar is detected when (gamma_1m)^2, (gamma_2m)^2, (gamma_3m)^2 are maximized to 1.

    Input:
        P, period of the pulsar
        a, semimajor axis of the pulsar
        e, eccentricity of the pulsar's orbit
        i, orbital inclination of the binary system
        T, duration of the observation
        . . .

    Returns:
        gamma_1m_sq, ratio of the highest power of the pulsar for the mth harmonic when acceleration and higher order terms are non-zero to when they are zero

        gamma_2m_sq, ratio of the highest power of the pulsar for the mth harmonic when jerk and higher order terms are non-zero to when they are zero

        gamma_3m_sq, ratio of the highest power of the pulsar for the mth harmonic when higher order terms than jerk are non-zero to when they are zero

    NOTE: Need to check units!!

    NOTE potential issues:
        - what should we use for m? (np.linspace(1, 20, 20) or something similar?)
        - what is omega_bar?

    """

    a_p_prime = ((((P/(2 * np.pi))**2) * const.G * (M_p + M_c))**(1/3)))*(M_c/(M_p+M_c))*np.sin(i)

    # get the mean anomaly
    M = w_0 * (t - T_p) # NOTE: what is w_0, T_p????
    M_0 = w_0 * T_p

    # get the eccentric anomaly
    E - e * np.sin(E) = M # NOTE: how do we want to solve this for E????
    E_0 - e * np.sin(E_0) = M_0

    # get the true anomaly at time t and time t=0
    f = 2 * np.arctan(np.sqrt((1 + e)/(1-e)) * np.tan(E/2))
    f_0 = 2 * np.arctan(np.sqrt((1 + e)/(1-e)) * np.tan(E_0/2))

    # get the radius vectors of the pulsar
    r_l = a_p_prime * (1 - e**2) * ((1 + e * np.cos(f))**(-1)) * np.sin(f + omega_bar)
    r_lo = a_p_prime * (1 - e**2) * ((1 + e * np.cos(f_0))**(-1)) * np.sin(f_0 + omega_bar)

    # get the velocity vector of the pulsar
    v_l = ((2 * np.pi)/P) * ((a_p_prime)/np.sqrt(1 - e**2)) * (np.cos(f + omega_bar) + (e * np.cos(omega_bar)))
    v_lo = ((2 * np.pi)/P) * ((a_p_prime)/np.sqrt(1 - e**2)) * (np.cos(f_0 + omega_bar) + (e * np.cos(omega_bar)))

    # get the acceleration vector of the pulsar
    a_l = - (((2 * np.pi)/P)**2) * (a_p_prime/((1 - e**2)**2)) * np.sin(f + omega_bar) * (1 + e * np.cos(f))**2
    a_lo = - (((2 * np.pi)/P)**2) * (a_p_prime/((1 - e**2)**2)) * np.sin(f_0 + omega_bar) * (1 + e * np.cos(f_0))**2

    # get the jerk vector of the pulsar
    j_l = - (((2 * np.pi)/P)**3) * (a_p_prime/((1 - e**2)**(7/2))) * ((1 + e*np.cos(f))**3) * (np.cos(f + omega_bar) + e * np.cos(omega_bar) - 3 * e * np.sin(f + omega_bar) * np.sin(f))
    j_l = - (((2 * np.pi)/P)**3) * (a_p_prime/((1 - e**2)**(7/2))) * ((1 + e*np.cos(f_0))**3) * (np.cos(f_0 + omega_bar) + e * np.cos(omega_bar) - 3 * e * np.sin(f_0 + omega_bar) * np.sin(f_0))


    # compute the gamma factors based on the above values for the pulsar
    gamma_1m = (1/(T * (-v_lo))) * np.abs(np.exp(((1j*m*w_p)/const.c)*(r_l - r_lo - (v_lo * T)))) # for the general case of orbital eccentricity, describes sensitivity loss of a standard pulsar search

    # Do we need these??
    # NOTE: still need to integrate the following two equaitons over T inside the abs!!!!
    gamma_2m = (1/(T)) * np.abs(np.exp(((1j*m*w_p)/const.c)*(r_l - r_lo - ((a_lo/2) * T**2) - (v_lo * T))))
    gamma_3m = (1/(T)) * np.abs(np.exp(((1j*m*w_p)/const.c)*(r_l - r_lo - ((j_lo/6) * T**3) - ((a_lo/2) * T**2) - (v_lo * T))))

    # return the gamma factors so we can check if the pulsar is detected with its given eccentricity
    # if the returned numbers (gamma_..**2) are 1, pulsar has been detected
    return gamma_1m**2, gamma_2m**2, gamma_3m**2

def S_min(M, P_orb, e, P, P_dot, B, x, y, z, vx, vy, vz, L, T_rec, d_f, n_chan, freq, tau_samp, G, t_int, npol = 2, SNmin = 10, beta=1):
    """
    Compute the minimum flux that a pulsar can have and still be detectable.

    Input:
        M, mass (M_sun)
        P_orb, orbital period (days)
        e, eccentricity
        P, rotational period (seconds)
        P_dot, change in rotational period (seconds)
        B, surface magnetic field (Tesla)
        x, y, z; cartesian coordinates
        vx, vy, vz; change in position of each cartesian
        L, luminosity (mJy Kpc^2)
        T_rec, receiver temperature
        d_f, receiver bandwidth (MHz)
        n_chan, number of channels in the survey
        freq, observing frequency (MHz)
        tau_samp, sampling time (seconds)
        G, gain of the telescope
        t_int, integration time (seconds)
        npol, number of polarizations in the detector (automatically set to 2)
        SNmin, minimum detection threshold (automatically set to 10)
        beta, parameter to account for the errors that increase the noise in the signal (automatically set to 1)

    Returns:
        S_min, the lower limit of flux a simulated pulsar can have to be detected at a given S/N ratio
        F, flux (mJy)
        D, distance (Kpc)
        gamma_1m_sq, gamma_2m_sq, gamma_3m_sq; numbers to check if the object is detectable with the given eccentricity
    """

    DM = DM_fnct(x, y, z) # dispersion measure in the direction of the pulsar, Units: pc/cm^3
    # correct the units to match those of NOTE DM_0 (m^-3??)
    DM = DM * const.parsec * ((1e2)**3) # new units, m^-3

    tau_scatt = tau_scatt_fnct(DM, freq) # ISM scattering time

    DM_0 = DM0_fnct(freq, d_f, n_chan, tau_samp) # diagonal dispersion measure of the survey

    T_sky = T_sky_fnct(x, y, z, freq) # get the sky temperature in the direction of the pulsar

    F, D = flux(L, x, y, z) # get the flux of the pulsar Units: F (mJy), D (Kpc)

    gamma_1m_sq, gamma_2m_sq, gamma_3m_sq = eccentricity(P, a, e, i, t_int, . . .)

    W_i = P*0.05 # fixed duty cycle in paper, seconds
    # compute the effective pulse width
    We = np.sqrt(W_i**2 + tau_samp**2 + (tau_samp*(DM/DM_0))**2 + tau_scatt**2) # Units: seconds

    # get the S/N ratio of the pulsar data
    npol = 2
    if We >= P:
        SNR = 0
    else: # NOTE the SNR ratio eq here does not match the north cap pulsar survey paper!!
        SNR = F / np.sqrt(np.pi / 2) / np.sqrt(We / (P - We)) / (T_rec + T_sky) * (G * np.sqrt(npol * d_f * t_int))

    # compute the minimum flux (S_min)
    S_min = beta*((SNmin*(T_rec+T_sky))/(G*np.sqrt(npol*t_int*(d_f/1e6))))*np.sqrt(We/(P-We))

    return S_min, F, D, gamma_1m, gamma_2m, gamma_3m

def f_beaming(p):
    """
    Computes the pulsar beaming fraction. It should return a result between 0 and 1.

    Input:
        p, period of the pulsar (seconds)
    Returns:
        f_b, pulsar beaming fraction
    """

    # beaming fraction model, P must be entered in seconds
    f_b = 0.09*np.log(P/10)**2+0.03 # is the log in this equation log base 10????

    if f_b>1 or f_b<0: # check to make sure that the beaming fraction is in the range that it should be within
        print("Error: beaming fraction is not within parameters. f_b =", f_b)
    return f_b
