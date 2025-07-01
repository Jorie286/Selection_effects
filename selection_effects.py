# radio selection effects complete script

# import necessary packages
import numpy as np
import scipy.constants as const
import galpy

# get functions and constants from other files
import gal_cart
import survey
import uconst

# allow the user to input the name of the survey, the luminosity, and the pulsar data file that we want to use
inputs = input("Enter input separated by spaces as follows: \nSurvey_Name Luminosity Pulsar_Data_File_Name: \n")

# take care of the error handling for the user inputs to make sure all are valid before we start computation
while True: # check to make sure enough inputs where given
    try:
        input_list = inputs.split() # inputs should be separated by spaces, split them into a list
        if len(input_list)!= 3:
            print("Error: not enough inputs given! Please try again. Expected 3 but got", len(input_list))
            inputs = input()
        break

# get the inputs in separate variables for easier error handling and use later
survey_name = input_list[0]
luinosity = input_list[1]
pulsar_data = input_list[2]

# do error checking on the inputs and get them corrected if necessary
while True:
    try: # check to make sure that the given survey name exsists in the surv_array
        s = survey.surv_array[np.where(survey.surv_array[0] == survey_name)] # get the row of data for the specific survey
        break
    except TypeError:
        print("Error: Survey name is not a string.")
        print("Use one of the following:")
        print(survey.surv_array[0])
        survey_name=input()
    except Exception as e: # try a different input if the first survey name did not work
        print("Error: Survey name does not match one that is avaiable.")
        print("Use one of the following:")
        print(survey.surv_array[0])
        survey_name=input()

while True:
    try:
        L = float(luminosity)
        if np.log(L)<-3 or np.log(L)>4:
            raise ValueError("Error: Luminosity is outside of accepted range. -3<= log(L) <= 4")
            print("Please try again.")
            luminosity = input()
        break
    except TypeError:
        print("Error: Invalid entry for luminosity, it must be of type float. Please try again.")
        luminosity=input()

while True:
    try:
        # NOTE: this will need to be changed depending on the type of file that we end up loading.
        p = np.fromfile(pulsar_data, dtype=np.float64) # load the pulsar data in as an array
        break
    except Exception as e:
        print("Error: invalid entry fro pulsar data file name, please try again.")
        pulsar_data = input()

print("All inputs successfuly processed!")

tsky1 = np.fromfile('tsky1.o', dtype=np.float64) # get the sky temperatures from the file


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

def S_min(p, s, L, npol = 2, beta=1):
    """
    Compute the minimum flux that a pulsar can have and still be detectable.

    Input:
        p, row of pulsar data from dataframe
        s, row number of the survey that we are using (full array is stored in survey.py)
        L, luminosity of the pulsar
        npol, number of polarizations in the detector (automatically set to 2)
        beta, parameter to account for the errors that increase the noise in the signal (automatically set to 1)

    Returns:
        S_min, the lower limit of flux a simulated pulsar can have to be detected at a given S/N ratio
        F, flux of the pulsar (W/kpc^2)
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
    S_min = beta*((SNR*(T_rec+T_sky))/(G*np.sqrt(npol*t_int*(d_f/1e6))))*np.sqrt(W_e/(P-W_e))
    return S_min, F

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
    L = . . . # luminosity of the pulsar
    F = L/(4*np.pi*(D**2)) # flux of the pulsar in (Watts / Kpc^2)
    return F

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


# NOTE: In what format do we want to save the data? The end of this file should be changed depending on the answer!!!
# create empty rows to store new data in (for dataframe ???)
pulsar_data["flux"] = None
pulsar_data["S_min"] = None
pulsar_data["is_detectable"] = None
pulsar_data["f_beaming"] = None

# iterate through each row of the simulated pulsar data and determine if the pulsar is detectable
for i, row in pulsar_data.iterrows:

    # compute S_min and flux
    S_min, flux = S_min(row, s, L)

    # get the galactic coordinates of the pulsar
    l, b, d = gal_cart.cart2gal(row['x'], row['y'], row['z'])

    if survey.s[-1](l, b)==1: # check to see if the pulsar is within the survey's viewing area. If it is, save the info.
        if flux >= S_min: # save info on whether or not the simulated pulsar would be detectable with the given survey as well as its flux and S_min
            pulsar_data["flux"][i] = flux
            pulsar_data["S_min"][i] = S_min
            pulsar_data["is_detectable"][i] = True
            pulsar_data["f_beaming"][i] = f_beaming(row)
        else:
            pulsar_data["flux"][i] = flux
            pulsar_data["S_min"][i] = S_min
            pulsar_data["is_detectable"][i] = False
            pulsar_data["f_beaming"][i] = f_beaming(row)
    else: # If the pulsar is not in the viewing area, don't save any info other than that it is not detectable
        pulsar_data["flux"][i] = None
        pulsar_data["S_min"][i] = None
        pulsar_data["is_detectable"][i] = False
        pulsar_data["f_beaming"][i] = None

# save the updated pulsar_data to a new csv file
pulsar_data.to_csv("pulsar_data_updated.csv", index=False)

