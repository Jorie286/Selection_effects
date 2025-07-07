# equivalent file of a binpuls____.c file

# import necessary packages
import numpy as np
import scipy.constants as const
import galpy

# get functions and constants from other files
import gal_cart
import survey
import uconst
import selection_effects

# allow the user to input the name of the survey, the luminosity, and the pulsar data file that we want to use
inputs = input("Enter input separated by spaces as follows: \nSurvey_Name Pulsar_Data_File_Name: \n")

# take care of the error handling for the user inputs to make sure all are valid before we start computation
while True: # check to make sure enough inputs where given
    input_list = inputs.split(" ") # inputs should be separated by spaces, split them into a list
    if len(input_list)!= 2:
        print("Error: not enough inputs given! Please try again. Expected 2 but got", len(input_list))
        inputs = input()
    else:
        break

# get the inputs in separate variables for easier error handling and use later
survey_name = input_list[0]
pulsar_data = input_list[1]

# do error checking on the inputs and get them corrected if necessary
while True:
    if str(survey_name) not in survey.surv_array[:, 0]:
        print("Error: Survey name does not match one that is avaiable.")
        print("Use one of the following:")
        print(survey.surv_array[:, 0])
        survey_name=input()
    else:
        if len(np.where(survey.surv_array[:, 0] == str(survey_name))[0])!=1:
            print("Which version of the survey would you like? Choose one of the following array indicies:")
            print(np.where(survey.surv_array[:, 0] == str(survey_name))[0])
            n = input()
            s = survey.surv_array[int(n)]
        else:
            s = survey.surv_array[np.where(survey.surv_array[:, 0] == str(survey_name))] # get the row of data for the specific survey
        break

while True:
    try:
        # NOTE: this will need to be changed depending on the type of file that we end up loading.
        p = np.fromfile(pulsar_data, dtype=np.float64) # load the pulsar data in as an array
        break
    except Exception as e:
        print("Error: invalid entry for pulsar data file name, please try again.")
        pulsar_data = input()

print("All inputs successfuly processed!")


# create empty rows to store new data in (for dataframe ???)
pulsar_data["flux"] = None
pulsar_data["S_min"] = None
pulsar_data["SNR"] = None
pulsar_data["f_beaming"] = None

# iterate through each row of the simulated pulsar data and determine if the pulsar is detectable
for i, row in pulsar_data.iterrows:

    # compute S_min and flux
    S_minDsq, flux, SNR, gamma_1m_sq, gamma_2m_sq, gamma_3m_sq = selection_effects.S_min(row, s, L)

    # get the galactic coordinates of the pulsar
    l, b, d = gal_cart.cart2gal(row['x'], row['y'], row['z'])

    if survey.s[-1](l, b)==1: # check to see if the pulsar is within the survey's viewing area. If it is, save the info.
        pulsar_data["flux"][i] = flux
        pulsar_data["S_min"][i] = S_minDsq
        pulsar_data["SNR"][i] = SNR
        pulsar_data["f_beaming"][i] = selection_effects.f_beaming(row)
        pulsar_data["gamma_1m"][i] = gamma_1m
        pulsar_data["gamma_2m"][i] = gamma_2m
        pulsar_data["gamma_3m"][i] = gamma_3m

# NOTE: In what format do we want to save the data?
# save the updated pulsar_data to a new csv file
pulsar_data.to_csv("out_file.csv", index=False)
