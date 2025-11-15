# Selection_effects

Determine which simulated pulsars would be detected in a survey by calculating the selection effects of the survey on the pulsar data.

To determine which pulsars would be detected, the binpulsBothFile.py should be run. It asks for inputs of what survey you want to use, a file containing the luminosities associated with the pulsar data, and the file in which the simulated pulsar data is stored. The input should be entered as follows:

$\verb|python binpulsBothFile.py survey_name input.dat system_type output.dat|$

Replace the survey_name with one listed in the avaliable surveys below. The input.dat file and output.dat file should be changed to the dat files that you want to use as the input and the one you want the output written to. The system_type should be changed to binary if you are using an input file that contains data from binary pulsar systems. Otherwise, it will assume you are using a single star system.

The code automatically impliments a soft cutoff of pulsar death lines and the radiometer equation to get the luminosity cutoff. It calculates the beaming fraction of the pulsar to get a better idea of the detectability of the pulsar. The pipeline impliments the correction for binary orbits given by Bagchi (Bagchi, M., et. al. 2013) based on the objects in the system and their orbital eccentricities.

The script will return the output file when it finishes running. This file will include the following new columns for each star or each star in the binary system:
- minimum observable luminosity (S_min1, S_min2)
- area (Area)
- sky temperature (T_sky)
- binary varible determining if the star was detected (without considering the beaming fraction) (det1, det2)
- beaming fraction (f_beaming1, f_beaming2)
- signal to noise ratio (SNR1, SNR2)

The units used within the code are:
- Dispersion Measure: pc cm^-3 (???)
- Distance: Kpc
- Frequency: MHz
- Flux: mJy
- Luminosity: mJy Kpc^2
- Orbital Period: days
- Rotation Period: seconds
- Surface Magnetic Field: Tesla
- Velocity: km/s

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
- SKA_2100 (t_int = 2100 s)
- MeerKat_GalPl
- MeerKat_tint
- MeerKat_Galpl_tint
- Lowlat
- TRUMP_Meer
- MeerKat
- MeerKat
- SKA_300 (t_int = 300 s)

The repositories and versions that are required to run this code include:
- astropy 6.1.3
- astropy-healpix 1.0.2
- galpy 1.7.2
- matplotlib 3.10.0
- numpy 1.26.4
- pandas 2.2.3
- python 3.10.18
- pygedm 3.3.0
- scipy 1.13.1
- skytempy (necessary files included in a folder within the repository)

