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
if len(sys.argv)!= 5:
    print("Error: not enough inputs given! Please try again. Expected 5 but got", len(sys.argv))
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

# create empty rows to store new data in
pulsar_data_out["S_min1*area"] = None
pulsar_data_out["S_min2*area"] = None
pulsar_data_out["Area"] = None
pulsar_data_out["T_sky"] = None
pulsar_data_out["f_beaming1"] = None
pulsar_data_out["f_beaming2"] = None
pulsar_data_out["det1"] = None
pulsar_data_out["det2"] = None
pulsar_data_out["SNR1"] = None
pulsar_data_out["SNR2"] = None
pulsar_data_out["alt_det1"] = None
pulsar_data_out["alt_det2"] = None

# iterate through each row of the simulated pulsar data and determine if the pulsar is detectable
for i, row in pulsar_data_out.iterrows():
    if file_type == "binary":

        # define each constant in the pulsar data as what it is for greater readability
        # Units: none, M_sun, M_sun, none, s, s, s^2, s^2, T, T, Kpc, Kpc, Kpc, Kpc, Kpc, Kpc, km/s, km/s, km/s, km/s, km/s, km/s, mJy Kpc^2, mJy Kpc^2
        ID, M_1, M_2, P_orb, e, a, P_1, P_2, P_dot1, P_dot2, B_1, B_2, x, y, z, vx, vy, vz, L1, L2, type1, type2 = row['ID'], row['m1(Msun)'], row['m2(Msun)'], row['porb(days)'], row['e'], row['a(AU)'], row['p1(s)'], row['p2(s)'], row['pdot1(s/s)'], row['pdot2(s/s)'], row['b1(T)'], row['b2(T)'], row['x(kpc)'], row['y(kpc)'], row['z(kpc)'], row['vx(km/s)'], row['vy(km/s)'], row['vz(km/s)'], row['l1(mJy kpc²)'], row['l2(mJy kpc²)'], row['type1'], row['type2']

        # check to make sure that the pulsar is not below the death lines, if it is spin period will be None
        det1 = selection_effects.death_lines(P_1, P_dot1, M_1, x, y, z, L1)
        det2 = selection_effects.death_lines(P_2, P_dot2, M_2, x, y, z, L2)

        # check to see if each object is a neutron star and if it is radio detectable, if it is not, save values as None, else compute the values
        if pd.isna(P_1) == True or det1==0:
            S_min1, flux1, area, f_b1, SNR1, T_sky = None, None, None, None, None, None
        else:
            S_min1, flux1, area, SNR1, T_sky = selection_effects.S_min(M_1, P_orb, e, a, P_1, P_dot1, B_1, x, y, z, vx, vy, vz, L1, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

            f_b1 = selection_effects.f_beaming(P_1)


        # repeat process for second star if the data includes binary systems
        if pd.isna(P_2) == True or det2==0:
            S_min2, flux2, area, f_b2, SNR2, T_sky = None, None, None, None, None, None
        else:
            S_min2, flux2, area, SNR2, T_sky = selection_effects.S_min(M_2, P_orb, e, a, P_2, P_dot2, B_2, x, y, z, vx, vy, vz, L2, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

            f_b2 = selection_effects.f_beaming(P_2)

        # get the Bagchi correction of the objects in the system
        alt_det1, alt_det2 = selection_effects.DNS_NSBH_sel_eff(P_1, P_2, P_orb, e, type1, type2)

    else:
        # define each constant in the pulsar data as what it is for greater readability
        # Units: none, M_sun, none, s, s^2, T, Kpc, Kpc, Kpc, km/s, km/s, km/s, mJy Kpc^2
        ID, M_1, P_orb, e, a, P_1, P_dot1, B_1, x, y, z, vx, vy, vz, L1 = row['ID'], row['m1(Msun)'], row['porb(days)'], row['e'], row['a'], row['p1(s)'], row['pdot1(s/s)'], row['b1(T)'], row['x(kpc)'], row['y(kpc)'], row['z(kpc)'], row['vx(km/s)'], row['vy(km/s)'], row['vz(km/s)'], row['l1(mJy kpc²)']

        # check to make sure that the pulsar is not below the death lines, if it is spin period will be None
        det1 = selection_effects.death_lines(P_1, P_dot1, x, y, z, L1)

        # check to see if each object is a neutron star, if it is not, save values as None, else compute the values
        if pd.isna(P_1) == True or det1==0:
            S_min1, flux1, area, f_b1, SNR1, T_sky = None, None, None, None, None, None
        else:
            S_min1, flux1, area, SNR1, T_sky = selection_effects.S_min(M_1, P_orb, e, a, P_1, P_dot1, B_1, x, y, z, vx, vy, vz, L1, T_rec, d_f, n_chan, freq, tau_samp, G, t_int)

            f_b1 = selection_effects.f_beaming(P_1)

        # set all other variables that would be returned for a binary to none for a single star system
        S_min2, flux2, gamma_1m_sq, gamma_2m_sq, gamma_3m_sq, f_b2, SNR2, alt_det1, alt_det2 = None, None, None, None, None, None, None, None, None


    # get the galacitc coordinates of the object
    l, b, d = gal_cart.cart2gal(x, y, z, degree=True)
    if s[-1](l, b) == 1 and pd.isna(S_min1)==False and pd.isna(S_min1)==False and pd.isna(area)==False and pd.isna(T_sky)==False and pd.isna(f_b1)==False and pd.isna(f_b2)==False and pd.isna(SNR1)==False and pd.isna(SNR2)==False: # check to see if the pulsar is within the survey's viewing area. If it is, save the info.
        pulsar_data_out.loc[i, "S_min1*area"] = S_min1*area
        pulsar_data_out.loc[i, "S_min2*area"] = S_min2*area
        pulsar_data_out.loc[i, "Area"] = area
        pulsar_data_out.loc[i, "T_sky"] = T_sky
        pulsar_data_out.loc[i, "f_beaming1"] = f_b1
        pulsar_data_out.loc[i, "f_beaming2"] = f_b2
        pulsar_data_out.loc[i, "det1"] = (L1 >= S_min1*area).astype(int)
        pulsar_data_out.loc[i, "det2"] = (L2 >= S_min2*area).astype(int)
        pulsar_data_out.loc[i, "SNR1"] = SNR1
        pulsar_data_out.loc[i, "SNR2"] = SNR2
        pulsar_data_out.loc[i, "alt_det1"] = alt_det1
        pulsar_data_out.loc[i, "alt_det2"] = alt_det2

# save the updated pulsar data to two new csv files
pulsar_data_out.to_csv(str(output_name), index=False, sep=";")
