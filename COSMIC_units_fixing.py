import pandas as pd
import numpy as np

import gal_cart

# read in the COSMIC output
table = pd.read_hdf("NSNS_sources.h5", "final_bpp")

# read in all the necessary table info
M_1, M_2, P_orb, e, a, omega_spin_1, omega_spin_2, B_1, B_2, x, y, z, vx, vy, vz, L1, L2, type1, type2 = table["mass_1"], table["mass_2"], table["porb"], table["ecc"], table["sep"], table["omega_spin_1"], table["omega_spin_2"], table["B_1"], table["B_2"], table["x"], table["y"], table["z"], table["v_x"], table["v_y"], table["v_z"], table["lum_1"], table["lum_2"], table["kstar_1"], table["kstar_2"]



# correcting units:
R_sun_to_AU = 6.9599e8 / 1.4960e11
a = a*R_sun_to_AU # a given in R_sun, multiply by (R_sun to m)/(m to AU)


# get the period from the angular velocity
yr_to_sec = 365.25*24*60*60 # conversion factor from years to seconds
P_1 = (2*np.pi)/(omega_spin_1*yr_to_sec) # convert angular velocity to period, seconds
P_2 = (2*np.pi)/(omega_spin_2*yr_to_sec)


# get pdot using the period and other parameters
R = 12*1e3 # m
mu_0 = 1.256e-6 # m kg s^-2 A^-2

I_1 = 0.237 * M_1 * (R**2) * 1.989e30* (1 + (4.2 * (M_1/R)) + 90*((M_1/R)**4))# kg m^2 (from Lorimer, D., et. al., Handbook of Pulsar Astronomy)
I_2 = 0.237 * M_2 * (R**2) * 1.989e30* (1 + (4.2 * (M_2/R)) + 90*((M_2/R)**4))# kg m^2 (from Lorimer, D., et. al., Handbook of Pulsar Astronomy)

alpha = 30*(np.pi/180)

omega_dot_1 = - (8*np.pi*(B_1**2)*(R**6)*((omega_spin_1 * yr_to_sec)**3)*(np.sin(alpha))**2) / (3*mu_0*(2.99e8**3)*I_1)
omega_dot_2 = - (8*np.pi*(B_2**2)*(R**6)*((omega_spin_2 * yr_to_sec)**3)*(np.sin(alpha))**2) / (3*mu_0*(2.99e8**3)*I_2)

Pdot_1 = - (omega_dot_1 * P_1)/(omega_spin_1 * yr_to_sec)
Pdot_2 = - (omega_dot_2 * P_2)/(omega_spin_2 * yr_to_sec)


# convert magnetic field in Gauss to Tesla
Gauss_to_Tesla = 1/1e4

B_1 = B_1*Gauss_to_Tesla
B_2 = B_2*Gauss_to_Tesla


# not sure about the units for x, y, z, vx, vy, vz; I would guess R_sun and R_sun/yr???
# put distance and velocity coordinates in to the proper units
R_sun_to_kpc = 6.96e8/(3.0857e16 * 1e3)
R_sunyr_to_kms = (6.96e8/1e3) / yr_to_sec

x, y, z = x*R_sun_to_kpc, y*R_sun_to_kpc, z*R_sun_to_kpc
vx, vy, vz = vx*R_sunyr_to_kms, vy*R_sunyr_to_kms, vz*R_sunyr_to_kms


# convert L_sun to mJy kpc^2
L_sun_to_mJykpc = []
for index in range(len(x)):
    l, b, d = gal_cart.cart2gal(x[index], y[index], z[index], degree=True) # cart2gal can't take a list or arrray so we need to iterate over the values.
    L_sun_to_mJykpc.append(3.826e33 / (7.4e27*(d**2)))

L1 = L1 * L_sun_to_mJykpc
L2 = L2 * L_sun_to_mJykpc

# make a list of null ID values to satisfy requirement in selection_effects.py
ID = np.zeros(np.shape(M_1))


# make a new dataframe to store the updated values that can be read directly into the selection effects code
new = {"ID":ID, "m1(Msun)":M_1, "m2(Msun)":M_2, "porb(days)":P_orb, "e":e, "a(AU)":a, "p1(s)":P_1, "p2(s)":P_2, "pdot1(s/s)":Pdot_1, "pdot2(s/s)":Pdot_2, "b1(T)":B_1, "b2(T)":B_2, "x(kpc)":x, "y(kpc)":y, "z(kpc)":z, "vx(km/s)":vx, "vy(km/s)":vy, "vz(km/s)":vz, "l1(mJy kpc²)":L1, "l2(mJy kpc²)":L2, "type1":type1, "type2":type2}

new_dataframe = pd.DataFrame(data = new)

new_dataframe.to_csv("COSMIC_converted.dat", index=False, sep=";")
