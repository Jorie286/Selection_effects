# radio selection effects functions

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
R0_Kpc=8.5 # Kpc, distance Sun-Galaxy Center

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

    l, b, d = gal_cart.cart2gal(x, y, z, degree=True) # Units: l (degrees), b (degrees), d (kpc)

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
    # NOTE: An equation for tau_scatt also appears in Pulsar Astronomy book pg 37 eq 3.9!!
    tau_scatt = -6.46 + 0.154 * np.log10(DM) + 1.07 * (np.log10(DM) ** 2) - 3.86 * np.log10(freq / 1e3)

    # tau_scatt from book:
    #tau_scatt=((DM/1000)**3.5)*((400/nu_MHz)**4)

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
    # printf("wavelength=%3.3f\n", % wavelength)

    DM_0 = 1000.0 * tau_samp *((3e2/wavelength)**3)/(8.3e6*(d_f/n_chan))
    # The following is dm0 copied from survey.c:
    # dm0 = 1000.0*s.t_samp*pow(3e2/wavelength,3)/(8.3e6*(s.receiverBW/s.n_chan ))
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

def E_fnct(E):
    """
    Define the equation to get the eccentric anomaly to be used in eccentricity() later.
    """
    return E - M - e * np.sin(E)

def eccentricity(P_orb, M_1, M_2, e, x, y, z, vx, vy, vz, i, T):
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
        G, gain of the telescope (Kelvin/ Jy) # NOTE: are these the units that gain is in in survey.py???
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
    """

    DM, tau_sc = DM_fnct(x, y, z, freq) # dispersion measure in the direction of the pulsar, Units: pc/cm^3, seconds
    #print("freq", freq)
    #print("DM", DM)

    tau_scatt = tau_scatt_fnct(DM, freq) # ISM scattering time, Units: seconds
    #print("tau_scatt", tau_scatt)

    DM_0 = DM0_fnct(freq, d_f, n_chan, tau_samp) # diagonal dispersion measure of the survey
    #print("DM_0", DM_0)

    #T_sky = T_sky_fnct(-x, -y, z, freq) # get the sky temperature in the direction of the pulsar (Kelvin)
    #print("T_sky", T_sky)
    T_sky=0 # set a temporary value for T_sky to see if there are other issues.

    F, D, Area = flux(L, x, y, z) # get the flux of the pulsar, Units: F (mJy), D (Kpc), Area (Kpc^2)
    #print("F", F)
    #print("Area", Area)

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
            # DEBATRI: is G in Kelvin/Jy???? If it is, then the calculations are good, else there is factor of 1e3 missing.
        S_min_list.append(S_min) # add each calculated S_min to the list

        # get the S/N ratio of the pulsar data
        npol = 2
        if We >= P:
            SNR = 0
        else: # DEBATRI the SNR ratio eq here does not match the north cap pulsar survey paper
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
