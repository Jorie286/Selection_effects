# Selection_effects

Determine which simulated pulsars would be detected in a survey by calculating the selection effects of the survey on the pulsar data.

To determine which pulsars would be detected, the binpulsBothFile.py should be run. It asks for inputs of what survey you want to use, a file containing the luminosities associated with the pulsar data, and the file in which the simulated pulsar data is stored. The script will return a pulsar_data_updated.csv file when it finishes running. This file will include new columns with it's flux, beaming fraction, flux lower limit for the given survey, and signal to noise ratio.

The units used within the code are:
- Angle: radians
- Dispersion Measure: m^-3 (???)
- Distance: Kpc
- Frequency: MHz
- Flux: mJy
- Luminosity: mJy Kpc^2
- Time/ Period: seconds

The surveys that are available to test for a pulsar's detectability include:
- Parkes_70_cm
- Parkes_Multi_Beam
- Swin_Interm Lat
- Swin_Extended
- Burgay_et_al
- GMRT
- Parkes_Multi_Beam_ALLSKY
- Parkes_Multi_Beam_part
- MeerKat
- SKA (t_int = 2100 s)
- MeerKat_GalPl
- MeerKat_tint
- MeerKat_Galpl_tint
- Lowlat
- TRUMP_Meer
- MeerKat
- MeerKat
- SKA (t_int = 300 s)

The repositories that are required to run this code include:
- galpy
- numpy
- pygedm
- scipy
- skytempy (necessary files included in a folder within the repository)
