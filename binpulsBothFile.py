# import necessary packages
import sys
import numpy as np
import pandas as pd
import scipy.constants as const
import galpy

# get functions and constants from other files
import gal_cart
import survey
import selection_effects

# take care of the error handling for the user inputs to make sure all are valid before we start computation
if len(sys.argv)!= 6:
    print("Error: not enough inputs given! Please try again. Expected 5 but got", len(sys.argv))
    print("Please pass, survey_name input_data file_type output_name")
    sys.exit(-1)

# get the inputs in separate variables for easier error handling and use later
survey_name = sys.argv[1]
pulsar_data = sys.argv[2]
file_type = sys.argv[3]
line_cutoff = sys.argv[4]
output_name = sys.argv[5]

# do error checking on the inputs and get them corrected if necessary
if str(survey_name) not in survey.surv_array[:, 0]:
    print("Error: Survey name does not match one that is available.")
    print("Use one of the following:")
    print(survey.surv_array[:, 0])
    sys.exit(-1)
else:
    s = survey.surv_array[np.where(survey.surv_array[:, 0] == str(survey_name))][0] # get the row of data for the specific survey

    # define each constant in the survey array as what it is for greater readability
    # units: none, Kelvin, MHz, none, MHz, s, none, s
    name, T_rec, d_f, n_chan, freq, tau_samp, G, t_int = s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]

try:
    p = pd.read_csv(pulsar_data, delimiter=';') # load the pulsar data in as a dataframe
except Exception as e:
    print("Error: invalid entry for pulsar data file name, please try again.")
    sys.exit(-1)

print("All inputs successfully processed!")

# duplicate the orignal dataframe so we don't modify the orignal
pulsar_data_out = p.copy()

# check to see if user wants to remove all data below death lines
if line_cutoff == "True": # if true remove points below death lines from the copied dataframe
    death1 = pulsar_data_out.loc[(10**(3.29 * np.log10(pulsar_data_out["p1(s)"]) - 16.55)) > pulsar_data_out["pdot1(s/s)"] or (10**(0.92 * np.log10(pulsar_data_out["p1(s)"]) - 18.65)) > pulsar_data_out["pdot1(s/s)"], "pdot1(s/s)"]

    death2 = pulsar_data_out.loc[(10**(3.29 * np.log10(pulsar_data_out["p2(s)"]) - 16.55)) > pulsar_data_out["pdot2(s/s)"] or (10**(0.92 * np.log10(pulsar_data_out["p2(s)"]) - 18.65)) > pulsar_data_out["pdot2(s/s)"], "pdot2(s/s)"]

    # remove spin period data that is below the death lines so that it is not included in the calculation and will not be plotted.
    pulsar_data_out.loc[death1, "p1(s)"] = None
    pulsar_data_out.loc[death2, "p2(s)"] = None


# create empty rows to store new data in
pulsar_data_out["S_min1_05*area"] = None
pulsar_data_out["S_min2_05*area"] = None
pulsar_data_out["S_min1_27*area"] = None
pulsar_data_out["S_min2_27*area"] = None
pulsar_data_out["S_min1_fwhm*area"] = None
pulsar_data_out["S_min2_fwhm*area"] = None
pulsar_data_out["Area"] = None
pulsar_data_out["T_sky"] = None
pulsar_data_out["f_beaming1"] = None
pulsar_data_out["f_beaming2"] = None
pulsar_data_out["det1_05"] = None
pulsar_data_out["det2_05"] = None
pulsar_data_out["det1_27"] = None
pulsar_data_out["det2_27"] = None
pulsar_data_out["det1_fwhm"] = None
pulsar_data_out["det2_fwhm"] = None
pulsar_data_out["SNR1_05"] = None
pulsar_data_out["SNR2_05"] = None
pulsar_data_out["SNR1_27"] = None
pulsar_data_out["SNR2_27"] = None
pulsar_data_out["SNR1_fwhm"] = None
pulsar_data_out["SNR2_fwhm"] = None
pulsar_data_out["alt_det1"] = None
pulsar_data_out["alt_det2"] = None

# iterate through each row of the simulated pulsar data and determine if the pulsar is detectable
for i, row in pulsar_data_out.iterrows():
    if file_type == "binary":

        # define each constant in the pulsar data as what it is for greater readability
        # Units: none, M_sun, M_sun, none, s, s, s^2, s^2, T, T, Kpc, Kpc, Kpc, Kpc, Kpc, Kpc, km/s, km/s, km/s, km/s, km/s, km/s, mJy Kpc^2, mJy Kpc^2
        ID, M_1, M_2, P_orb, e, a, P_1, P_2, P_dot1, P_dot2, B_1, B_2, x, y, z, vx, vy, vz, L1, L2 = row['ID'], row['m1(Msun)'], row['m2(Msun)'], row['porb(days)'], row['e'], row['a(AU)'], row['p1(s)'], row['p2(s)'], row['pdot1(s/s)'], row['pdot2(s/s)'], row['b1(T)'], row['b2(T)'], row['x(kpc)'], row['y(kpc)'], row['z(kpc)'], row['vx(km/s)'], row['vy(km/s)'], row['vz(km/s)'], row['l1(mJy kpc²)'], row['l2(mJy kpc²)']

        # check to see if each object is a neutron star, if it is not, save values as None, else compute the values
        if pd.isna(P_1) == True:
            S_min1_05, S_min1_27, S_min1_fwhm, flux1, area, f_b1, SNR1_05, SNR1_27, SNR1_fwhm, T_sky = None, None, None, None, None, None, None, None, None, None
        else:
            S_min1_05, S_min1_27, S_min1_fwhm, flux1, area, SNR1_05, SNR1_27, SNR1_fwhm, T_sky = selection_effects.S_min(M_1, P_orb, e, a, P_1, P_dot1, B_1, x, y, z, vx, vy, vz, L1, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

            f_b1 = selection_effects.f_beaming(P_1)

        # repeat process for second star if the data includes binary systems
        if pd.isna(P_2) == True:
            S_min2_05, S_min2_27, S_min2_fwhm, flux2, area, f_b2, SNR2_05, SNR2_27, SNR2_fwhm, T_sky = None, None, None, None, None, None, None, None, None, None
        else:
            S_min2_05, S_min2_27, S_min2_fwhm, flux2, area, SNR2_05, SNR2_27, SNR2_fwhm, T_sky = selection_effects.S_min(M_2, P_orb, e, a, P_2, P_dot2, B_2, x, y, z, vx, vy, vz, L2, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

            f_b2 = selection_effects.f_beaming(P_2)

        #alt_det1, alt_det2 = selection_effects.DNS_NSBH_sel_eff(P1, P2, P_orb1, P_orb2, e, sys_type)

    else:
        # define each constant in the pulsar data as what it is for greater readability
        # Units: none, M_sun, none, s, s^2, T, Kpc, Kpc, Kpc, km/s, km/s, km/s, mJy Kpc^2
        ID, M_1, P_orb, e, a, P_1, P_dot1, B_1, x, y, z, vx, vy, vz, L1 = row['ID'], row['m1(Msun)'], row['porb(days)'], row['e'], row['a'], row['p1(s)'], row['pdot1(s/s)'], row['b1(T)'], row['x(kpc)'], row['y(kpc)'], row['z(kpc)'], row['vx(km/s)'], row['vy(km/s)'], row['vz(km/s)'], row['l1(mJy kpc²)']

        # check to see if each object is a neutron star, if it is not, save values as None, else compute the values
        if pd.isna(P_1) == True:
            S_min1_05, S_min1_27, S_min1_fwhm, flux1, area, f_b1, SNR1_05, SNR1_27, SNR1_fwhm, T_sky = None, None, None, None, None, None, None, None, None, None
        else:
            S_min1_05, S_min1_27, S_min1_fwhm, flux1, area, SNR1_05, SNR1_27, SNR1_fwhm, T_sky = selection_effects.S_min(M_1, P_orb, e, a, P_1, P_dot1, B_1, x, y, z, vx, vy, vz, L1, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

            f_b1 = selection_effects.f_beaming(P_1)

        # set all other variables that would be returned for a binary to none for a single star system
        S_min2_05, S_min2_27, S_min2_fwhm, flux2, L2, gamma_1m_sq, gamma_2m_sq, gamma_3m_sq, f_b2, SNR2_05, SNR2_27, SNR2_fwhm, alt_det1, alt_det2 = None, None, None, None, None, None, None, None, None, None, None, None, None, None


    # get the galacitc coordinates of the object
    l, b, d = gal_cart.cart2gal(x, y, z, degree=True)
    if s[-1](l, b) == 1: # check to see if the pulsar is within the survey's viewing area. If it is, save the info.
        pulsar_data_out.loc[i, "S_min1_05*area"] = S_min1_05*area
        pulsar_data_out.loc[i, "S_min2_05*area"] = S_min2_05*area
        pulsar_data_out.loc[i, "S_min1_27*area"] = S_min1_27*area
        pulsar_data_out.loc[i, "S_min2_27*area"] = S_min2_27*area
        pulsar_data_out.loc[i, "S_min1_fwhm*area"] = S_min1_fwhm*area
        pulsar_data_out.loc[i, "S_min2_fwhm*area"] = S_min2_fwhm*area
        pulsar_data_out.loc[i, "Area"] = area
        pulsar_data_out.loc[i, "T_sky"] = T_sky
        pulsar_data_out.loc[i, "f_beaming1"] = f_b1
        pulsar_data_out.loc[i, "f_beaming2"] = f_b2
        pulsar_data_out.loc[i, "det1_05"] = (L1 >= S_min1_05*area).astype(int)
        pulsar_data_out.loc[i, "det2_05"] = (L2 >= S_min2_05*area).astype(int)
        pulsar_data_out.loc[i, "det1_27"] = (L1 >= S_min1_27*area).astype(int)
        pulsar_data_out.loc[i, "det2_27"] = (L2 >= S_min2_27*area).astype(int)
        pulsar_data_out.loc[i, "det1_fwhm"] = (L1 >= S_min1_fwhm*area).astype(int)
        pulsar_data_out.loc[i, "det2_fwhm"] = (L2 >= S_min2_fwhm*area).astype(int)
        pulsar_data_out.loc[i, "SNR1_05"] = SNR1_05
        pulsar_data_out.loc[i, "SNR2_05"] = SNR2_05
        pulsar_data_out.loc[i, "SNR1_27"] = SNR1_27
        pulsar_data_out.loc[i, "SNR2_27"] = SNR2_27
        pulsar_data_out.loc[i, "SNR1_fwhm"] = SNR1_fwhm
        pulsar_data_out.loc[i, "SNR2_fwhm"] = SNR2_fwhm
        pulsar_data_out.loc[i, "alt_det1"] = alt_det1
        pulsar_data_out.loc[i, "alt_det2"] = alt_det2

# save the updated pulsar data to two new csv files
pulsar_data_out.to_csv(str(output_name), index=False, sep=";")
