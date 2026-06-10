interact=0
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import unique,vstack,Table,Column
import paths
import re 
groupnight=1
import sys
import os
import sys
import argparse

parser = argparse.ArgumentParser(                    
	prog='clean_asassn_lightcurve',
                    description='Takes raw downloaded ASASSN lightcurves and cleans the photometry',
                    epilog='Use -d to see the intermediate plots and analysis',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-o", "--outdir", help="Directory to write output file",default='.')
parser.add_argument("-n", "--name", help="name of the object",default="TEST")
parser.add_argument("input",help="Lightcurve file")
parser.add_argument("-d","--display", help="plot lightcurves",
                    action="store_true",default=False)
args = parser.parse_args()

if args.display:
	print("Plotting ON")
	interact=1

if args.outdir:
    print("Output directory:", args.outdir)

fin = args.input
print(f'Input file is {fin}')


if os.access(args.outdir, os.W_OK):
	print(f'{args.outdir} is writeable')
else:
	print(f'{args.outdir} is NOT writeable. Stopping.')
	quit()

filename_suffix='ecsv'
fout = os.path.join(args.outdir, args.name + '_ASASSN_clean.' + filename_suffix)

print(f'Output file is {fout}')

mindelt = 0.015 # days minimum between photometric epochs to call it a distinct separate observation

# read in the ASASSN light curve exported as CSV from the web interface
try:
	ta = ascii.read(fin,format='csv',guess=False,
		converters={
		'HJD': float,
		'Camera': str,
		'FWHM': float,
		'Limit': float,
		'mag': str,
		'mag_err': float,
		'flux(mJy)': float,
		'flux_err': float,
		'Filter': str})
except FileNotFoundError:
    print(f"Error: The file {fin} does not exist.")

ta.remove_column('UT Date')
ta.add_column(1,name='npoints') 

# 'mag' is read in as a str initially because it can have ">1234.234" representing as an upper limit
# we remove these in the goodmag expression below, copy the table and reassign the dtype to make it float

#      HJD           UT Date       Camera FWHM Limit   mag   mag_err flux(mJy) flux_err Filter
# ------------- ------------------ ------ ---- ------ ------ ------- --------- -------- ------
# 2457420.65322 2016-02-02.1500246     be 1.46 17.458  13.45   0.005    15.995     0.08      V
# bad values are labelled 99.99 in mag,mag_err, flux,flux_err

print(f'# table {fin} read in with {len(ta)} lines')

goodmag = [False if re.match('>',l) else True for l in ta['mag']]
print(f'# {np.sum(goodmag)} lines with good magnitude values')

t = ta[goodmag]
t['mag'] = Column(t['mag'],dtype=float)

mbad = t['mag']>99

print(f'# {np.sum(mbad)} lines with bad (>99) magnitude values')

t=t[~mbad]

print(f'# {len(t)} lines remaining')

t['MJD'] = t['HJD']-2400000.5

if interact:
	fig, (ax) = plt.subplots(1,1,figsize=(12,6))
	ax.errorbar(t['MJD'],t['flux(mJy)'],yerr=t['flux_err'],fmt='.')
	ax.set_ylabel('Flux [mJy]')
	ax.set_xlabel('Epoch [MJD]')
	ax.set_title('data from {}'.format(fin))
	#plt.show()

t_unique = unique(t,keys='Filter')

if len(t_unique) > 1:
	fig, axes = plt.subplots(1,len(t_unique),figsize=(12,6))
	ax=np.ndarray.flatten(axes)
else:
	fig, axes = plt.subplots(1,1,figsize=(12,6))
	axes = np.array([axes])

if len(t_unique) > 1:
	fig2, axes2 = plt.subplots(1,len(t_unique),figsize=(12,6))
	ax2=np.ndarray.flatten(axes2)
else:
	fig2, axes2 = plt.subplots(1,1,figsize=(12,6))
	axes2 = np.array([axes2])

# make a blank copy of the input table
tout = Table(t,copy=True)
tout.remove_rows(np.arange(len(tout))) # there MUST be a better way to initialize a blank Table using another Table column and dtype...

def mean_cluster(t):
	# take a Table with ASASSN data and return a single table with the mean of the time, 
	to = Table(t,copy=True)
	to.remove_rows(np.arange(len(to)-1)) # remove all but one row
	to['MJD'][0] = np.mean(t['MJD'])
	to['HJD'][0] = np.mean(t['HJD'])
	to['FWHM'][0] = np.mean(t['FWHM'])

	to['mag'][0] = np.mean(t['mag'])
	to['flux(mJy)'][0] = np.mean(t['flux(mJy)'])

	# def ms(f,e): # testing with Lyons 1991 data analysis for physical science studentse
	# 	e2 = e*e 
	# 	meanf = np.sum(f/e2)/np.sum(1/e2)
	# 	meane = 1. / np.sqrt(np.sum(1/e2))
	# 	return meanf,meane
	# print(ms(t['mag'].value,t['mag_err'].value))
	# print(np.std(t['flux(mJy)'].value)/np.sqrt(len(t)))
	# print(ms(t['flux(mJy)'].value,t['flux_err'].value))

	# The error bar is the larger of either (i) the STD of the two or three points or (ii) the mean of the ASASSN quoted errors. 
	# This is to avoid anomalously low error bars if the two points happen to be very close to each other in value but have very large reported error bars.

	mag_err_std = np.std(t['mag'])
	mag_err_asas = np.mean(t['mag_err'])

	if mag_err_std > mag_err_asas:
		to['mag_err'][0] = mag_err_std
	else:
		to['mag_err'][0] = mag_err_asas

	flux_err_std = np.std(t['flux(mJy)'])
	flux_err_asas = np.mean(t['flux_err'])

	if flux_err_std > flux_err_asas:
		to['flux_err'][0] = flux_err_std
	else:
		to['flux_err'][0] = flux_err_asas

	return to

# split by filters, if there's more than one
for (filt,ax,ax2) in zip(t_unique['Filter'],axes,axes2):
	print(f'# working on filter {filt}')
	mfilt = (t['Filter']==filt)
	tfilt = t[mfilt] 
	print(f'# with {len(tfilt)} rows')

	# additionally split by ASASSN camera name:
	tcamera_unique = unique(tfilt,keys='Camera')

	for cam in tcamera_unique['Camera']:
		mcamera = (tfilt['Camera']==cam)
		tfc = tfilt[mcamera]

		if interact:
			ax.errorbar(tfc['MJD'],tfc['flux(mJy)'],yerr=tfc['flux_err'],fmt='.',label=cam)
			ax.set_ylabel('Flux [mJy]')
			ax.set_xlabel('Epoch [MJD]')
			ax.set_title(f'data from filter {filt}')


		if groupnight: # make averages over clusters of data points
			is_sorted = lambda a: np.all(a[:-1] <= a[1:])
			# first check and sort table by increasing HJD
#			print(is_sorted(tfc['MJD']))
			deltat = np.diff(tfc['MJD'])
			if interact:
				ax2.hist(deltat,bins=201,range=(0,1),log=True, label=cam)
				ax2.set_ylabel('N')
				ax2.set_xlabel('delta time [days]')
				ax2.set_title(f'data from filter {filt}')


			k = (deltat>mindelt) # True for large delta times
			# prepend a zero, so that you have a True when you've just had a big jump in time
			k = np.concatenate(([0],k)) # prepend a 0 on this
			# cumulative sum this, so that 0, 1, 2, ... represent each a single cluster of observations that should be combined together.
			kcumul = np.cumsum(k)

			# crazy.... I want to know how many points per cluster of observations, should be mostly 3 or 2, occasionally 1. It's okay if you
			# get the occasional cluster of 6 points in there!
			r = np.histogram(kcumul,bins=np.max(kcumul)+1,range=(-0.5,np.max(kcumul)+1)) # this is insane... how many points within a given range?
			print(f'# camera {cam}')
			for cluster in np.arange(np.max(kcumul)+1):
				t_cluster = tfc[(kcumul==cluster)]
				if len(t_cluster) > 1:
					teff = mean_cluster(t_cluster)
					teff['npoints'][0]=len(t_cluster)
					tout.add_row(teff[0])
				else:
					t_cluster['npoints'][0]=1
					tout.add_row(t_cluster[0])

	if interact:
		ax.errorbar(tout['MJD'],tout['flux(mJy)'],yerr=tout['flux_err'],color='orange',fmt='s',alpha=0.5,elinewidth=5)
		ax.legend()
		ax2.legend()


if interact:

	if len(t_unique) > 1:
		fig3, axes3 = plt.subplots(1,len(t_unique),figsize=(12,6))
		ax3=np.ndarray.flatten(axes3)
	else:
		fig3, axes3 = plt.subplots(1,1,figsize=(12,6))
		axes3 = np.array([axes3])

	# split by filters, if there's more than one
	for (filt,ax3) in zip(t_unique['Filter'],axes3):
		print(f'# working on filter {filt}')
		mfilt = (tout['Filter']==filt)
		tfilt = tout[mfilt] 
		print(f'# with {len(tfilt)} rows')

		# additionally split by ASASSN camera name:
		tcamera_unique = unique(tfilt,keys='Camera')

		ax3.errorbar(tfilt['MJD'],tfilt['flux(mJy)'],yerr=tfilt['flux_err'],fmt='.')
		ax3.set_ylabel('Flux [mJy]')
		ax3.set_xlabel('Epoch [MJD]')
		ax3.set_title(f'data from filter {filt}')

if interact:
	plt.show()

# this convoluted way of writing out an ECSV table is because the ecsv writer in astropy 
# does not respect the .format values in the Table object, but plain ascii does...
# I write out a plain ascii table with the correct formatting, read that in, and write that out as an ECSV table.

tout['HJD'].format = '%14.6f'
tout['FWHM'].format = '%4.2f'
tout['Limit'].format = '%6.3f'
tout['mag'].format = '%6.3f'
tout['mag_err'].format = '%6.4f'
tout['flux(mJy)'].format = '%8.3f'
tout['flux_err'].format = '%6.2f'
tout['MJD'].format = '%12.6f'

tout.write('tmp.txt', format='ascii',overwrite=True)

tb = ascii.read('tmp.txt',format='basic',guess=False,	converters={
	'HJD': float,
	'Camera': str,
	'FWHM': float,
	'Limit': float,
	'mag': float,
	'mag_err': float,
	'flux(mJy)': float,
	'flux_err': float,
	'Filter': str,
	'npoints': int})

tb.write(fout, format='ecsv', overwrite=True)
#tb.write(sys.stdout, format='ecsv', overwrite=True)
print(f'# {fout} written out cleaned ASASSN photometry table')
