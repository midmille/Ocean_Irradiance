
# need to import this (as a minimum) from the
# OOI python directory (make sure path correct)
from ooi_data_explorations.data_request import data_request

# variables needed for download
# following are from OOI group
site="CE02SHSP"
assembly="profiler"
method="recovered_cspp"
# OOI gang suggest April - August 2021
# should have valid data
start="2021-04-01"
stop="2021-04-30"
# not necessary, but the relevatnt deployments (per OOI gang)
# are 19-23
#deploy=19

#there are four instruments to give data
instrument="OPTAA"	# OPTAA spectrophotometer
# https://oceanobservatories.org/instrument-class/optaa/
			# - contains wavelength-dependent info.
			# - beam_attenuation
			#   - may need to calculate coefficient manually?
			# - optical_absorption
			# - wavelength_a & _c
			# - temperature & salinity
instrument="SPKIR"	# spectral irradiance
# https://oceanobservatories.org/instrument-class/spkir/
			# - spkir_abj_cspp_downwelling_vector	
			# - vin_sense? va_sense?
			# - depth
			# - profiler_timestamp / internal_timestamp
instrument="FLORT"	# 3-stream fluorometer
# https://oceanobservatories.org/instrument-class/fluor/
			# - chl-a fluorescence (fluorometric_chlorophyll_a)
 			# - optical backscatter (red wavelengths)
			# - CDOM
			# - temperature & salinity
			# - depth
			# - profiler_timestamp / internal_timestamp
instrument="PARAD"	# photosynthetically active radiation (PAR) 
# https://oceanobservatories.org/instrument-class/parad/
			# - parad_j_par_counts_output
			# - depth
			# - profiler_timestamp / internal_timestamp

# TODO: find if there is CTD information??
# TODO: howto access the buoy data??


# how to download data
data = data_request(site, assembly, instrument, method, start=start, stop=stop)

# display all of the variable names
data.data_vars.keys()

# display specific variable info
# including units, etc. 
data.variables['NAME']

# get variable data
data.variables['NAME'].data

# Probably good to follow Jupyter example for condensing
# arry into only useful variables...
# https://nbviewer.org/github/oceanobservatories/ooi-data-explorations/blob/master/python/examples/notebooks/optaa/process_ooinet_optaa.ipynb
