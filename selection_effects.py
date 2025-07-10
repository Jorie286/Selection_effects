# radio selection effects functions

# import necessary packages
import numpy as np
import scipy.constants as const
from scipy.optimize import fsolve
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

    # d needs to be in pc for the pygedm function
    DM, tau_sc = pygedm.dist_to_dm((l, b), dist = d, nu=(freq/1e3)) # Units: DM (pc / cm^3), tau_sc (GHz)
    tau_sc = tau_sc*1e3 # convert tau_sc to MHz timescale

    # DEBATRI need to correct the units to match those of DM_0 (s^2 m^-3?)
    DM = DM * const.parsec * ((1e2)**3) # m^-3

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
    # DEBATRI what units are used in this equation? It seems to be a scaling relation for tau_scatt so it probably requires specific units for it to work properly but they are not stated where it is defined as t_scatter_fnct2 in pulsar_survey_functions.c.
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

    # DEBATRI the units of DM_0 look like s^2 m^-3, should they be m^-3?
    DM_0 = 1000.0 *tau_samp*((3e2/wavelength)**3)/(8.3e6*(d_f/n_chan))

    # The following is dm0 copied from survey.c:
    # dm0 = 1000.0*s.t_samp*pow(3e2/wavelength,3)/(8.3e6*(s.receiverBW/s.n_chan ))
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
    # DEBATRI what are the units of Tsky?
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
    F = L/(4*np.pi*(D**2)) # flux of the pulsar in (mJy)
    return F, D

def E_fnct(E):
    """
    Define the equation to get the eccentric anomaly to be used in eccentricity() later.
    """
    return E - M - e * np.sin(E)

def eccentricity(P_orb, M_1, M_2, e, x1, y1, z1, x2, y2, z2, vx1, vy1, vz1, vx2, vy2, vz2, i, T):
    """
    Compute the effects eccentricity in a binary orbit has on the detectability of a pulsar. A pulsar is detected when (gamma_1m)^2, (gamma_2m)^2, (gamma_3m)^2 are maximized to 1.

    Input:
        P_orb, orbital period of the binary (days)
        M_1, mass of first object in binary (M_sun)
        M_2, mass of second object in binary (M_sun)
        e, eccentricity of the object's orbit
        x1, y1, z1; cartesian coordinates of the first object (Kpc)
        x2, y2, z2; cartesian coordinates of the second object (Kpc)
        vx1, vy1, vz1; velocity of first object in cartesian coordinates (km/s)
        vx2, vy2, vz2; velocity of second object in cartesian coordinates (km/s)
        i, orbital inclination of the binary system
        T, duration of the observation (seconds)
        . . .

    Returns:
        gamma_1m_sq, ratio of the highest power of the pulsar for the mth harmonic when acceleration and higher order terms are non-zero to when they are zero

        gamma_2m_sq, ratio of the highest power of the pulsar for the mth harmonic when jerk and higher order terms are non-zero to when they are zero

        gamma_3m_sq, ratio of the highest power of the pulsar for the mth harmonic when higher order terms than jerk are non-zero to when they are zero

    DEBATRI potential issues (marked by NOTE):
        - line 179: do we need to consider the orbital inclination of the system
        - line 184, 185: how do we get the epoch of periastron passage (T_p)?
        - lines 184, 216: what should we use for m and t? (np.linspace(1, 20, 20) or something similar?)
        - lines 190, 191: will E_fnct work with M not being defined within its equation?
        - line 198: how to get the longitude of the periastron?
        - line 216: what should be used for T, the duration of the observation
        - lines 220, 221: would the efficiency and relative efficiency be useful? they will be difficult to get due to integrals.

    """
    # get the distance between the two objects in the binary orbit
    r = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) # Units: Kpc

    # get the relative velocity of the two objects in the binary orbit
    v = np.sqrt((vx1-vx2)**2 + (vy1-vy2)**2 + (vz1-vz2)**2) # Units: km/s

    # define the gravitational constant in more convenient units for this function
    G_units = 4.301e-9 * 1e3 # Units: km^2 Kpc M_sun^-1 s^-2

    # get the semimajor axis of the binary system with the vis-visa equation
    a = -(((v**2)/(G_units * (M_1 + M_2))) - (2/r))**(-1) # Units: Kpc

    # NOTE: do we need to consider the orbital inclination of the system?
    a_p_prime = ((((P/(2 * np.pi))**2) * G_units * (M_1 + M_2))**(1/3))*(M_2/(M_1 + M_2))*np.sin(i) # Units: km^(2/3) Kpc^(1/3)

    w_0 = (2*np.pi)/P # get the orbital angular frequency (Units: Hz)

    # get the mean anomaly
    M = w_0 * (t - T_p) # NOTE: T_p is the epoch of periastron passage??? what to do with the time variable???
    M_0 = w_0 * T_p

    # solve for the eccentric anomaly numerically
    # NOTE: I am not sure if fsolve will run the following without issue because M is not defined in E_funct
    E_guess = 1 # may need to be changed depending on how large or small E is going to be
    E = fsolve(E_fnct, E_guess)
    E_0 = fsolve(E_fnct, E_guess)

    # get the true anomaly at time t and time t=0
    f = 2 * np.arctan(np.sqrt((1 + e)/(1-e)) * np.tan(E/2))
    f_0 = 2 * np.arctan(np.sqrt((1 + e)/(1-e)) * np.tan(E_0/2))

    # get the radius vectors of the pulsar
    r_l = a_p_prime * (1 - e**2) * ((1 + e * np.cos(f))**(-1)) * np.sin(f + omega_bar) # NOTE: omega_bar is the longitude of the periastron???
    r_lo = a_p_prime * (1 - e**2) * ((1 + e * np.cos(f_0))**(-1)) * np.sin(f_0 + omega_bar)

    # get the velocity vector of the pulsar
    v_l = ((2 * np.pi)/P) * ((a_p_prime)/np.sqrt(1 - e**2)) * (np.cos(f + omega_bar) + (e * np.cos(omega_bar)))
    v_lo = ((2 * np.pi)/P) * ((a_p_prime)/np.sqrt(1 - e**2)) * (np.cos(f_0 + omega_bar) + (e * np.cos(omega_bar)))

    # get the acceleration vector of the pulsar
    a_l = - (((2 * np.pi)/P)**2) * (a_p_prime/((1 - e**2)**2)) * np.sin(f + omega_bar) * (1 + e * np.cos(f))**2
    a_lo = - (((2 * np.pi)/P)**2) * (a_p_prime/((1 - e**2)**2)) * np.sin(f_0 + omega_bar) * (1 + e * np.cos(f_0))**2

    # get the jerk vector of the pulsar
    j_l = - (((2 * np.pi)/P)**3) * (a_p_prime/((1 - e**2)**(7/2))) * ((1 + e*np.cos(f))**3) * (np.cos(f + omega_bar) + e * np.cos(omega_bar) - 3 * e * np.sin(f + omega_bar) * np.sin(f))
    j_lo = - (((2 * np.pi)/P)**3) * (a_p_prime/((1 - e**2)**(7/2))) * ((1 + e*np.cos(f_0))**3) * (np.cos(f_0 + omega_bar) + e * np.cos(omega_bar) - 3 * e * np.sin(f_0 + omega_bar) * np.sin(f_0))


    # compute the gamma factors based on the above values for the pulsar
    # NOTE: what to do about T (duration of the observation)??? Same as tau_samp or t_int???
    gamma_1m = (1/(T * (-v_lo))) * np.abs(np.exp(((1j*m*w_p)/const.c)*(r_l - r_lo - (v_lo * T)))) # for the general case of orbital eccentricity, describes sensitivity loss of a standard pulsar search

    # Do we need these??
    # NOTE: still need to integrate the following two equaitons over T inside the abs!!!!
    gamma_2m = (1/(T)) * np.abs(np.exp(((1j*m*w_p)/const.c)*(r_l - r_lo - ((a_lo/2) * T**2) - (v_lo * T)))) # efficiency factor
    gamma_3m = (1/(T)) * np.abs(np.exp(((1j*m*w_p)/const.c)*(r_l - r_lo - ((j_lo/6) * T**3) - ((a_lo/2) * T**2) - (v_lo * T)))) # relative efficiency factor

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

    tau_scatt = tau_scatt_fnct(DM, freq) # ISM scattering time

    DM_0 = DM0_fnct(freq, d_f, n_chan, tau_samp) # diagonal dispersion measure of the survey

    T_sky = T_sky_fnct(x, y, z, freq) # get the sky temperature in the direction of the pulsar

    F, D = flux(L, x, y, z) # get the flux of the pulsar (Units: F (mJy), D (Kpc))

    W_i = P*0.05 # fixed duty cycle in paper (Units: seconds)
    # compute the effective pulse width
    We = np.sqrt(W_i**2 + tau_samp**2 + (tau_samp*(DM/DM_0))**2 + tau_scatt**2) # Units: seconds

    # get the S/N ratio of the pulsar data
    npol = 2
    if We >= P:
        SNR = 0
    else: # DEBATRI the SNR ratio eq here does not match the north cap pulsar survey paper!!
        SNR = F / np.sqrt(np.pi / 2) / np.sqrt(We / (P - We)) / (T_rec + T_sky) * (G * np.sqrt(npol * d_f * t_int))

    # compute the minimum flux (S_min)
    # DEBATRI what are the units of S_min??? Flux is in mJy, S_min will not match this if we multiply by D**2
    S_min = beta*((SNmin*(T_rec+T_sky))/(G*np.sqrt(npol*t_int*(d_f/1e6))))*np.sqrt(We/(P-We))

    return S_min, F, D, gamma_1m, gamma_2m, gamma_3m

def f_beaming(P):
    """
    Computes the pulsar beaming fraction. It should return a result between 0 and 1.

    Input:
        P, period of the pulsar (seconds)
    Returns:
        f_b, pulsar beaming fraction
    """

    # beaming fraction model, P must be entered in seconds
    f_b = 0.09*np.log10(P/10)**2+0.03

    return f_b
