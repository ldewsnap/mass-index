import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
# Suppress a warning in filter_data() when giving coords as a list instead of a
# wavelet file.
pd.options.mode.chained_assignment = None

# From https://github.com/wmpg/WesternMeteorPyLib
import wmpl
from wmpl.Utils.TrajConversions import *


# Radar station coordinates:
LAT = 43.264690*np.pi/180.0
LON = -80.771940*np.pi/180.0

# Define directories within the project folder.
multinest_dir = 'Multinest/'  # Path to Multinest
wavelet_dir = 'Data/Wavelet_files/'  # Path to wavelet files
datafile_dir = 'Output/Multinest_datafiles/'  # Path to datafiles to feed into Multinest
fit_dir = 'Output/Multinest_fits/'  # Path to Multinest fit output
radarData_dir = 'Output/CMOR_data/'  # Path to local radar data storage
plot_dir = 'Output/Figures/'  # Path to output of plotting commands
CMOR_dir = '/mnt/meteor/radar/spool/pt0-29/'  # Where the server was mounted
local_pt0_dir = 'Data/pt0_files/'  # path to locally-saved pt0 files, if any


### Function definitions ###

def pull_CMOR_data(year, slon, local=False):
	""" Pull CMOR data for a given year and solar longitude from the WMPG server.
	Inputs:
		year -- integer
		slon -- solar longitude, integer. File will contain observations from slon
			to slon+1
	Optional input:
		local -- Boolean. If False, pulls from server mounted to CMOR_dir. If True,
			looks for file in local_pt0_dir. Default is False.
	Returns a pandas dataframe containing the data. If the file is not found,
		returns None.

	example: df = pull_CMOR_data(2019,76)
	"""
	path = local_pt0_dir if local else CMOR_dir
	filename = path+'mev-'+str(year)+'-'+str(slon).zfill(3)+'-29-00.pt0'
	
	try: 
		# pt0 file headers are inconsistent and I've seen at least three variations. Hopefully this block will catch all of them.
		with open(filename, 'r') as f:
			for line in f:
				if line.startswith('#'):
					header = line
				else:
					break
		header = header[2:].strip().split()
		header=[sub.replace('flags', 'fl') for sub in header]
		df = pd.read_csv(filename, comment='#', delim_whitespace=True, header=None, names=header)
	except (IOError, FileNotFoundError) as e:
		print('No pt0 file found for solar longitude {} in {}'.format(str(slon).zfill(3), year))
		return None
	
	# Convert separate date and time columns into a single datetime column
	df['datetime'] = pd.to_datetime(df.date.map(str) + ' ' + df.time)
	df['JulianDate'] = df['datetime'].apply(datetime2JD)
	
	# The file stores missing values as strings of periods, of inconsistent length.
	# Replace them with NaN, and convert the columns that had them.
	df.replace(to_replace=r'\.{2,}', value=np.nan, regex=True, inplace=True)
	df = df.astype({'new_ptn':'float64', 'd_new_ptn+':'float64',
		'd_new_ptn-':'float64', 'wind':'float64', 'd_wind':'float64'})
	return df


def filter_data(df, coords, vel=False, vel_threshold=0.15, radiant_threshold=5):
	""" Filter for range, data quality, and radiant, and optionally velocity.
	The range filter is applied per Blaauw et al (2011)
	doi:10.1111/j.1365-2966.2010.18038.x

	Required inputs:
		df -- pandas dataframe containing unfiltered data.
		coords -- string or list. If a string, should be name of shower wavelet
			file. If a list, should be radiant coordinates in the form (ll0, beta).
	Optional inputs:
		vel -- either a float or False. Gives mean shower velocity in km/s for
			filtering, or if False, disables velocity filtering. Default is False.
		vel_threshold -- float. Gives velocity threshold for filtering, if used.
		Default is 0.15 (i.e. 15%).
		radiant_threshold -- float. Gives angular separation threshold in degrees
		to make a shower association. Default is 5 degrees.

	Returns a pandas dataframe with only the data that passed the filtering.

	Example call:
	filter_data(df, 'Wavelet-ZPE.dat', vel=27.5, vel_threshold=0.2, radiant_threshold=4)
	filter_data(df, (344.4, 4.5))
	"""
	df = df[(df.range > 110) & (df.range < 130) 
		& (df.th < 70)  # Unphysical, given the range cut.
		& (df.fl == 0)]  # Bad interferometry if fl=1

	# Velocity cut
	if vel:
		df = df[(df.new_ptn > vel*(1-vel_threshold)) 
		& (df.new_ptn < vel*(1+vel_threshold))]

	# Identify the shower radiant
	if isinstance(coords, str):
		wavelet = read_wavelet(coords)
		df = df.apply(get_wavelet_radiant, axis=1, args=(wavelet,))
	else:
		df['radiant_ll0'] = coords[0]
		df['radiant_beta'] = coords[1]

	df['separation'] = df.apply(check_radiant, axis=1)
	df.drop(['radiant_ll0', 'radiant_beta'], axis=1, inplace=True)
	df_shower = df[df['separation'] <= radiant_threshold]
	return df_shower


def check_radiant(row):
	""" Compare the position of a meteor event to a given shower radiant.
	For a meteor to be compatible with a certain radiant, it must have been
	observed in a direction perpendicular to the direction of the radiant.
	Inputs:
		row -- a row of a pandas dataframe from the other functions.
	Returns:
		separation -- float, gives angular separation of meteor event from a
			plane normal to the given radiant.
	"""
	# Convert sun-centred ecliptic longitude to standard ecliptic
	lam = np.mod(row['radiant_ll0']+row['solar'], 360)

	[radiant_RA, radiant_dec] = ecliptic2RaDec(row['JulianDate'],
		lam*np.pi/180.0, row['radiant_beta']*np.pi/180.0)
	# Convert CMOR's format to AltAz
	azimuth = (90-row['phi'])*np.pi/180.0
	altitude = (90-row['th'])*np.pi/180.0
	[event_RA, event_dec] = altAz2RADec(azimuth, altitude, row['JulianDate'],
		LAT, LON)
	# Find angular separation between radiant and meteor event
	separation = np.abs(np.abs(np.arccos(
		np.cos(np.pi/2.0 - radiant_dec)*np.cos(np.pi/2.0 - event_dec) 
		+ np.sin(np.pi/2.0 - radiant_dec)*np.sin(np.pi/2.0 - event_dec
			)*np.cos(radiant_RA-event_RA)
		))*180.0/np.pi - 90)

	return separation


def read_wavelet(filename):
	"""Read in a wavelet file as a pandas dataframe.
	Input:
		filename -- String; name of wavelet file.
	Returns:
		wavelet_df -- pandas dataframe containing wavelet data.
	"""
	wavelet_df = pd.read_csv(wavelet_dir+filename, comment='#', index_col=None,
		delim_whitespace=True, header=None, names=['solar','ll0','beta','xsig'], usecols=[0,1,2,7])
	wavelet_df.drop(wavelet_df.tail(3).index,inplace=True)
	return wavelet_df


def get_wavelet_radiant(row, wavelet_df):
	"""Get the radiant at a given solar longitude, to account for drift.
	If the solar longitude is outside range of file, uses the radiant at the
	peak of shower activity. This is done to allow for comparisons to sporadics at
	other times of year.
	Inputs:
		row -- a row of a pandas dataframe from the other functions.
		wavelet_df -- wavelet dataframe for shower
	Returns:
		row -- original row with added columns containing the coordinates of the
			shower radiant at the time.
	"""
	wavelet_row = wavelet_df[wavelet_df.solar == row['solar']//1]
	if len(wavelet_row.index):
		row['radiant_ll0'] = wavelet_row.ll0.values[0]
		row['radiant_beta'] = wavelet_row.beta.values[0]
		return row
	else:
		peak_row = wavelet_df.iloc[wavelet_df.xsig.idxmax()]
		row['radiant_ll0'] = peak_row.ll0
		row['radiant_beta'] = peak_row.beta
		return row


def save_datafile(df_shower, shower_name):
	""" Save a shower dataframe as a datafile suitable for Multinest.
	Input:
		df_shower -- dataframe containing shower data.
		shower_name -- string; name or abbreviation for shower. Must not
			contain spaces.
	"""
	with open(datafile_dir+'datafile_'+shower_name, 'w') as file:
		for amplitude in sorted(df_shower['maxamp']):
			file.write("%s\n" % amplitude)


def run_multinest(shower_name):
	"""Run Multinest fit on the datafile for the selected shower.
	Input:
		shower_name -- string; name or abbreviation for shower. Must not
			contain spaces.
	"""
	os.system('cp '+datafile_dir+'datafile_'+shower_name+' '+multinest_dir
		+'datafile')

	# Multinest's mass.sh fails unless actually run from its directory
	working_dir = os.getcwd()
	os.chdir(multinest_dir)
	os.system('sh mass.sh')
	os.chdir(working_dir)

	os.rename(multinest_dir + 'gr.out', fit_dir + 'fit_' + shower_name + '.out')


def plot_fit(shower_name, axis=plt.gca(), color='k', linestyle='-',
	shower_fullname=False):
	""" Plot the amplitudes and fit for a given data set.
	Input:
		shower_name -- string; name or abbreviation for shower. Must not contain
			spaces.
	Optional inputs:
		axis -- pyplot axis object on which to plot. Defaults to current axis.
		color -- color to plot the shower. Default is 'k'
		linestyle -- for plotting fit. Default is '-'
		shower_fullname -- String or False. Sets a separate name for the shower in
			the legend compared to filenames (e.g. full name versus abbreviation).
	Example:
	plot_fit('ZPE', color='g', linestyle='--', shower_fullname='Zeta Perseids')
	"""
	if not shower_fullname: shower_fullname=shower_name
	amplitudes, counts = read_datafile(shower_name)
	slope, slope_err_plus, slope_err_minus, y_int = read_multinest_file(shower_name)
	axis.scatter(amplitudes, counts, c=color, s=1)
	# Plot power law fit from Multinest. 
	points = (3e2, 5e4) # Draw line a bit past the turnoff points
	axis.plot(points, [(10**y_int)*(x**(1-slope)) for x in points],
		color=color, linestyle=linestyle,
		label=shower_fullname+', s= %.2f +%.2f/-%.2f' % (slope, slope_err_plus,
			slope_err_minus))


def read_datafile(shower_name):
	"""Load Multinest input datafile for given data set.
	Input:
		shower_name -- string; name or abbreviation for shower. Must not
			contain spaces.
	Returns:
		amplitudes - list of floats; amplitude values from the datafile.
		count -- list of ints; reverse ordered cumulative count of amplitudes.
	"""
	with open(datafile_dir+'datafile_'+shower_name) as file:
		amplitudes = [float(line) for line in file]
	count = [x+1 for x in range(len(amplitudes))[::-1]]
	return amplitudes, count


def read_multinest_file(shower_name):
	""" Load Multinest fit output file for a given data set.
	Input:
		shower_name -- string; name or abbreviation for shower. Must not
			contain spaces.
	Returns:
		slope -- slope of fit (i.e. mass index)
		slope_err -- error opf the slope
		y_int -- intercept of fit (physically irrelevant, but used for plots)
	"""
	with open(fit_dir+'fit_'+shower_name+'.out') as file:
		line = file.readline().split()
		slope = 1-float(line[0])
		slope_err_plus = -float(line[1])
		slope_err_minus = float(line[2])
		y_int = float(file.readline().split()[0])
	return slope, slope_err_plus, slope_err_minus, y_int


def process_shower(shower_name, years, solar_longitudes, coords, vel=False,
	vel_threshold=0.15, radiant_threshold=5, local=False):
	"""Process the data, and produce a filtered data file and a Multinest fit.
	This function will loop over the given solar longitudes and years, pulling the
	data either from the server or a local folder, and filtering it according to the
	input parameters, then combines each loop's filtered data into a single
	dataframe. This dataframe is saved and a Multinest-compatible datafile created,
	then Multinest is run to produce a mass index fit.
	Inputs:
		shower_name -- string; name or abbreviation for shower. Must not
			contain spaces.
		years -- list; years for which to get data.
		solar_longitudes -- list; solar longitudes for which to get data.
		coords -- string or list. If a string, should be name of shower wavelet
			file. If a list, should be radiant coordinates in the form (ll0, beta).
	Optional inputs:
		vel -- either a float or False. Gives mean shower velocity in km/s for
			filtering, or if False, disables velocity filtering. Default is False.
		vel_threshold -- float. Gives velocity threshold for filtering, if used.
		Default is 0.15 (i.e. 15%).
		radiant_threshold -- float. Gives angular separation threshold in degrees
		for radiant filtering. Default is 5 degrees.
		local -- Boolean. If False, pulls from server. If True, looks for file in
			local_pt0_dir. Default is False.
	Returns:
		df_shower_full -- pandas dataframe; the filtered shower data.
	"""
	df_list = []  # Stores dataframes for each part of the loop
	for y in years:
		for s in solar_longitudes:
			df = pull_CMOR_data(y, s, local=local)
			if df is not None:
				df_shower_segment = filter_data(df, coords, vel=vel,
					vel_threshold=vel_threshold,
					radiant_threshold=radiant_threshold)
				df_list.append(df_shower_segment)
				print('Done {}, year {}, solar longitude {}'
					.format(shower_name, y,str(s).zfill(3)))
	df_shower_full = pd.concat(df_list)
	df_shower_full.to_csv(radarData_dir+'shower_data_'+shower_name, mode='a',
		header=True, index=False)
	save_datafile(df_shower_full, shower_name)
	run_multinest(shower_name)
	return df_shower_full


### Testing and example usage ###
# Below are a few examples used to test the functions, which demonstrate usage.

test_mode = 2
if test_mode in [1,2]:
	# Compare a shower to the sporadic activity at the same radiant, during a
	# period when the shower is not active, and fit both populations, then plot
	# the two on the same figure to compare.
	if test_mode == 1:
		# Short test with no velocity filtering and without a wavelet file.
		years = [2019]
		vel = False
		coords = [207.5, 11]
	elif test_mode == 2:
		# Test with several years, with a velocity cut, using the wavelet file, and
		# with some data deliberately missing.
		years = range(2013,2020)
		vel = 34.2
		coords ='Wavelet-GEM.dat'
	process_shower('GEM', years, range(254,266), coords=coords, vel=vel, local=True)
	print('Done shower. Starting sporadics.')
	process_shower('GEM_sporadic', years, range(230,240), coords=coords, vel=vel,
		local=True)
	fig, ax = plt.subplots(1,1,figsize=(8,6))
	ax.set_yscale('log')
	ax.set_xscale('log')
	plot_fit('GEM', axis=ax, color='k', linestyle='-', shower_fullname='Geminids')
	plot_fit('GEM_sporadic', axis=ax, color='grey', linestyle=':',
		shower_fullname='Sporadic')
	ax.legend(loc='lower left')
	ax.set_ylim(1, ax.get_ylim()[1])
	ax.set_ylabel('Cumulative count')
	ax.set_xlabel('Amplitude (digital units)')
	plt.tight_layout()
	plt.savefig(plot_dir+'test'+str(test_mode)+'.png')

elif test_mode == 3:
	# Fit mass index at each bin of 1 degree solar longitude, and plot the change
	# over time. Should produce a visible dip as the shower reaches its peak.
	# Takes several minutes.
	slopes = []
	errs_plus = []
	errs_minus = []
	slons_range = range(254,266)
	for slon in slons_range:
		process_shower('GEM_'+str(slon), range(2013,2020), [slon],
			coords='Wavelet-GEM.dat', vel=34.2, local=True)
		fit = read_multinest_file('GEM_'+str(slon))
		slopes.append(fit[0])
		errs_plus.append(fit[1])
		errs_minus.append(fit[2])
	fig, ax = plt.subplots(1,1,figsize=(7,4))
	ax.errorbar(range(254,266), slopes, yerr=(errs_minus, errs_plus), marker='o',
		markersize=9, capsize=5)
	ax.set_xlabel('Solar longitude (deg)')
	ax.set_ylabel('Mass index')
	plt.tight_layout()
	plt.savefig(plot_dir+'test'+str(test_mode)+'.png')

elif test_mode == 4:
	# Simply demonstrates that this function works without local files.
	print(pull_CMOR_data(2019, 1, local=False).head())