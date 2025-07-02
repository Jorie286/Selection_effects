# radio selection effects functions

# import necessary packages
import numpy as np
import scipy.constants as const

# get functions and constants from other files
import gal_cart
import survey

# now that any errors should be fixed, define functions that we will need to determine the selection effects
def DM_fnct(x, y, z):
    """
    Get the dispersion measure of the pulsar.

    Inputs:
    -------
    x, y, z; cartisian coordinates for the position of the pulsar

    Returns:
    --------
    DM, the dispersion measure of hte pulsar in the given direction.
    """

    l, b, d = gal_cart.cart2gal(x, y, z)
    neg_one = -1
    DM, limit, sm, smtau, smtheta = dmdsm(l, b, neg_one, d)
    return DM

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
    tau_scatt = -6.46 + 0.154 * np.log10(DM) + 1.07 * np.log10(DM) ** 2 - 3.86 * np.log10(freq / 1e3)
    tau_scatt = 10**tau_scatt # ms
    tau_scatt = tau_scatt * 1e-3 # convert to seconds
    return tau_scatt

def DM0_fnct(freq, d_f, n_chan, tau_samp):
    """
    Get the diaganal dispersion measure of the survey.

    Input:
    ------
    freq, survey frequency (seconds)
    d_f, reciever bandwidth (MHz)
    n_chan, number of reciever chanels
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
    x, y, x; cartisian coordinates of the pulsar
    freq, survey frequency (MHz)

    Returns:
    --------
    T_sky, sky temperature in the direction of the puslar
    """
    wavelength = const.c / (freq * 1e6) # m
    l, b, d = gal_cart.cart2gal(x, y, z) # galactic coordinates of the puslar
    tsky400 = tsky1(l, b)  # Sky temperature at 408 MHz NOTE: this line of the function may need to be fixed depending on how the tsky1.o file is formated.
    lambda_408 = const.c / (408 * 1e6)
    T_sky = tsky400 * (wavelength / lambda_408) ** 2.8
    return T_sky

def flux(L, x, y, z):
    """
    Compute the flux of the pulsar.

    Input:
        L, luminosity of the modeled pulsar, an input given by user (Watts ???)
        p, row of pulsar data from dataframe

    Returns:
        F, flux of the pulsar (W/Kpc^2)
    """

    D = np.sqrt(x**2 + y**2 + z**2) # distance to pulsar (Kpc)
    F = L/(4*np.pi*(D**2)) # flux of the pulsar in (Watts / Kpc^2)
    return F

def S_min(p, s, L, npol = 2, SNmin = 10, beta=1):
    """
    Compute the minimum flux that a pulsar can have and still be detectable.

    Input:
        p, row of pulsar data from dataframe
        s, row number of the survey that we are using (full array is stored in survey.py)
        L, luminosity of the pulsar
        npol, number of polarizations in the detector (automatically set to 2)
        SNmin, minimum detection threshold (automatically set to 10)
        beta, parameter to account for the errors that increase the noise in the signal (automatically set to 1)

    Returns:
        S_min, the lower limit of flux a simulated pulsar can have to be detected at a given S/N ratio
        F, flux of the pulsar (W/kpc^2)
        SNR, signal to noise ratio of the pulsar
    """

    # define each constant in the survey array as what it is for greater readability
    # units: none, ???, MHz, none, MHz, s, none, s
    name, T_rec, d_f, n_chan, freq, tau_samp, G, t_int = s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]

    # define each constant in the pulsar data as what it is for greater readability
    P, x, y, z = p['P'], p['x'], p['y'], p['z']# P in seconds, x, y, z in Kpc

    DM = DM_fnct(x, y, z)# dispersion measure in the direction of the pulsar (how to compute???)

    tau_scatt = tau_scatt_fnct(DM, freq) # ISM scattering time

    DM_0 = DM0_fnct(freq, d_f, n_chan, tau_samp) # diagonal dispersion measure of the survey

    T_sky = T_sky_fnct(x, y, z, freq) # get the sky temperature in the direction of the pulsar

    F = flux(L, x, y, z) # get the flux of the pulsar

    W_i = P*0.05 # fixed duty cycle in paper????, seconds
    # compute the effective pulse width
    W_e=np.sqrt(W_i**2 + tau_samp**2 + (tau_samp*(DM/DM_0))**2 + tau_scatt**2)

    # get the S/N ratio of the pulsar data
    npol = 2
    if We >= P:
        SNR = 0
    else:
        SNR = F / np.sqrt(np.pi / 2) / np.sqrt(We / (P - We)) / (T_rec + T_sky) * (G * np.sqrt(npol * d_f * t_int))

    # compute the minimum flux
    S_min = beta*((SNmin*(T_rec+T_sky))/(G*np.sqrt(npol*t_int*(d_f/1e6))))*np.sqrt(W_e/(P-W_e))
    return S_min, F, SNR

def f_beaming(p):
    """
    Computes the pulsar beaming fraction. It should return a result between 0 and 1.

    Input:
        p, row of the pulsar data from dataframe (or array???)
    Returns:
        f_b, pulsar beaming fraction
    """
    P = p['P'] # seconds

    # beaming fraction model, P must be entered in seconds
    f_b = 0.09*np.log(P/10)**2+0.03 # is the log in this equation log base 10????

    if f_b>1 or f_b<0: # check to make sure that the beaming fraction is in the range that it should be within
        print("Error: beaming fraction is not within parameters. f_b =", f_b)
    return f_b
