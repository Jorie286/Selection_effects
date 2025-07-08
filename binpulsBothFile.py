# equivalent file of a binpuls____.c file

# import necessary packages
import sys
import numpy as np
import pandas as pd
import scipy.constants as const
import galpy

# get functions and constants from other files
import gal_cart
import survey
import uconst
import selection_effects

# take care of the error handling for the user inputs to make sure all are valid before we start computation
if len(sys.argv)!= 5:
    print("Error: not enough inputs given! Please try again. Expected 4 but got", len(sys.argv))
    print("Please pass, survey_name input_data file_type output_name")
    sys.exit(-1)

# get the inputs in separate variables for easier error handling and use later
survey_name = sys.argv[1]
pulsar_data = sys.argv[2]
file_type = sys.argv[3]
output_name = sys.argv[4]

# do error checking on the inputs and get them corrected if necessary
if str(survey_name) not in survey.surv_array[:, 0]:
    print("Error: Survey name does not match one that is available.")
    print("Use one of the following:")
    print(survey.surv_array[:, 0])
    sys.exit(-1)
else:
    s = survey.surv_array[np.where(survey.surv_array[:, 0] == str(survey_name))] # get the row of data for the specific survey

    # define each constant in the survey array as what it is for greater readability
    # units: none, ???, MHz, none, MHz, s, none, s
    name, T_rec, d_f, n_chan, freq, tau_samp, G, t_int = s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]

try:
    p = pd.read_csv(pulsar_data) # load the pulsar data in as a dataframe
except Exception as e:
    print("Error: invalid entry for pulsar data file name, please try again.")
    sys.exit(-1)

print("All inputs successfully processed!")


# create empty rows to store new data in
pulsar_data["flux1"] = None
pulsar_data["flux2"] = None
pulsar_data["S_min1"] = None
pulsar_data["S_min2"] = None
pulsar_data["f_beaming1"] = None
pulsar_data["f_beaming2"] = None
pulsar_data["gamma_1m_sq1"] = None
pulsar_data["gamma_2m_sq1"] = None
pulsar_data["gamma_3m_sq1"] = None
pulsar_data["gamma_1m_sq2"] = None
pulsar_data["gamma_2m_sq2"] = None
pulsar_data["gamma_3m_sq2"] = None
pulsar_data["d_1"] = None
pulsar_data["d_2"] = None

# iterate through each row of the simulated pulsar data and determine if the pulsar is detectable
for i, row in p.iterrows:
    if file_type == "binary":

        # define each constant in the pulsar data as what it is for greater readability
        # Units: none, M_sun, M_sun, none, s, s, s^2, s^2, T, T, Kpc, Kpc, Kpc, Kpc, Kpc, Kpc, km/s, km/s, km/s, km/s, km/s, km/s, mJy Kpc^2, mJy Kpc^2
        ID, M_1, M_2, P_orb, e, P_1, P_2, P_dot1, P_dot2, B_1, B_2, x1, y1, z1, x2, y2, z2, vx1, vy1, vz1, vx2, vy2, vz2, L1, L2 = row['ID'], row['M1'], row['M2'], row['Porb'], row['e'], row['P1'], row['P2'], row['Pdot1'], row['Pdot2'], row['B1'], row['B2'], row['x1'], row['y1'], row['z1'], row['x2'], row['y2'], row['z2'], row['vx1'], row['vy1'], row['vz1'], row['vx2'], row['vy2'], row['vz2'], row['L1'], row['L2']

        # check to see if each object is a neutron star, if it is not, save values as None, else compute the values
        if P_1 == "NaN":
            S_min1, flux1, d_1, gamma_1m_sq1, gamma_2m_sq1, gamma_3m_sq1 = None, None, None, None, None, None
        else:
            S_min1, flux1, d_1, gamma_1m_sq1, gamma_2m_sq1, gamma_3m_sq1 = selection_effects.S_min(M_1, P_orb, e, P_1, P_dot1, B_1, x1, y1, z1, vx1, vy1, vz1, L1, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

        if P_2 == "NaN":
            S_min2, flux2, d_2, gamma_1m_sq2, gamma_1m_sq2, gamma_3m_sq2 = None, None, None, None, None, None
        else:
            S_min2, flux2, d_2, gamma_1m_sq2, gamma_2m_sq2, gamma_3m_sq2 = selection_effects.S_min(M_2, P_orb, e, P_2, P_dot2, B_2, x2, y2, z2, vx2, vy2, vz2, L2, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

    else:
        # define each constant in the pulsar data as what it is for greater readability
        # Units: none, M_sun, none, s, s^2, T, Kpc, Kpc, Kpc, km/s, km/s, km/s, mJy Kpc^2
        ID, M_1, P_orb, e, P_1, P_dot1, B_1, x1, y1, z1, vx1, vy1, vz1, L1 = row['ID'], row['M1'], row['Porb'], row['e'], row['P1'], row['Pdot1'], row['B1'], row['x1'], row['y1'], row['z1'], row['vx1'], row['vy1'], row['vz1'], row['L1']

        # check to see if each object is a neutron star, if it is not, save values as None, else compute the values
        if P_1 == "NaN":
            S_min1, flux1, d_1, gamma_1m_sq1, gamma_2m_sq1, gamma_3m_sq1 = None, None, None, None, None, None
        else:
            S_min1, flux1, d_1, gamma_1m_sq1, gamma_2m_sq1, gamma_3m_sq1 = selection_effects.S_min(M_1, P_orb, e, P_1, P_dot1, B_1, x1, y1, z1, vx1, vy1, vz1, L1, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

        # set all other variables that would be returned for a binary to none for a single star system
        S_min2, flux2, gamma_1m_sq2, gamma_2m_sq2, gamma_3m_sq2, d_2 = None, None, None, None, None, None


    if survey.s[-1](l, b) == 1: # check to see if the pulsar is within the survey's viewing area. If it is, save the info.
        pulsar_data["flux1"][i] = flux1
        pulsar_data["flux2"][i] = flux2
        pulsar_data["S_min1"][i] = S_min1
        pulsar_data["S_min2"][i] = S_min2
        pulsar_data["f_beaming1"][i] = selection_effects.f_beaming(P_1)
        pulsar_data["f_beaming2"][i] = selection_effects.f_beaming(P_2)
        pulsar_data["gamma_1m_sq1"][i] = gamma_1m_sq1
        pulsar_data["gamma_2m_sq1"][i] = gamma_2m_sq1
        pulsar_data["gamma_3m_sq1"][i] = gamma_3m_sq1
        pulsar_data["gamma_1m_sq2"][i] = gamma_1m_sq2
        pulsar_data["gamma_2m_sq2"][i] = gamma_2m_sq2
        pulsar_data["gamma_3m_sq2"][i] = gamma_3m_sq2
        pulsar_data["d_1"][i] = d_1
        pulsar_data["d_2"][i] = d_2

# save the updated pulsar_data to a new csv file
pulsar_data.to_csv(output_name, index=False)
