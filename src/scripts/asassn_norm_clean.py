import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import unique,vstack,Table,Column
import paths
import re 
from makphot import *

import os
import sys
import argparse

parser = argparse.ArgumentParser(                    
	prog='asassn_norm_clean',
                    description='takes cleaned ASASSN light curve and adds columns for normalised and cleaned photometry',
                    epilog='Use -d to see the intermediate plots and analysis',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#parser.add_argument("-n", "--name", help="name of the object",default="TEST")
parser.add_argument("input",help="Cleaned ASASSN lightcurve file")
parser.add_argument("-d","--display", help="plot lightcurves",
                    action="store_true",default=True)
args = parser.parse_args()

fin = args.input
print(f'# Input file is {fin}')

if args.display:
	print("# Plotting ON")
	interact=1

t = ascii.read(fin)
print(t.info)

# get a list of the unique bandpasses
t_by_filter = t.group_by('Filter')
t.add_column('fnorm')
t.add_column('fnorm_err')
print('all observed photometric bands:')
print(t_by_filter.groups.keys)

fig, axx = plt.subplots(1,1,figsize=(12,6))
for filters, tf in zip(t_by_filter.groups.keys, t_by_filter.groups):
	filt = filters['Filter']
	axx.errorbar(tf['MJD'],tf['flux(mJy)'],yerr=tf['flux_err'],label=filt,fmt='.')

print(t_by_filter)

if interact:
	plt.show()

quit()


fig, axes = plt.subplots(len(t_by_filter.groups.keys),3,figsize=(12,6))
ax = np.ndarray.flatten(axes)

for ind, (filters, tf) in enumerate(zip(t_by_filter.groups.keys, t_by_filter.groups)):
	filt = filters['Filter']
	axes[ind][0].hist(tf['flux_err'],bins=51,label=filt)
	axes[ind][0].legend()

	(cum_req, x, cum) = hist_errors (tf['flux_err'])

	sig_mask = tf['flux_err']<(cum[sig_percentile])
	tfs=tf[sig_mask]
	axx.errorbar(tfs['MJD'],tfs['flux(mJy)'],yerr=tfs['flux_err'],fmt='s')

	axes[ind][1].plot(cum,x)


plt.show()

quit()


wefig, (ax) = plt.subplots(1,1,figsize=(12,6))



#tminV, tmaxV = 58150,58450

#(V_flux_norm, V_flux_norm_err)= \
#    mean_rms_region(tc_binV[mV], fc_binV[mV], binned_V_err[mV], tminV, tmaxV)

tming, tmaxg = 58150,58450

(g_flux_norm, g_flux_norm_err)= \
    mean_rms_region(tc_bing[mg], fc_bing[mg], binned_g_err[mg], tming, tmaxg)

#print(f'Flux V band normalised {V_flux_norm:5.2}\pm{V_flux_norm_err:5.2f}')
print(f'Flux g band normalised {g_flux_norm:5.2}\\pm{g_flux_norm_err:5.2f}')

#tVout = Table([tc_binV[mV],fc_binV[mV]/V_flux_norm,
#    binned_V_err[mV]/V_flux_norm],
#    names=('MJD','fnorm','fnormerr'))
#tVout['Filter'] = 'V'

tgout = Table([tc_bing[mg],fc_bing[mg]/g_flux_norm,
    binned_g_err[mg]/g_flux_norm],
    names=('MJD','fnorm','fnormerr'))
tgout['Filter'] = 'g'

fig, (ax) = plt.subplots(1,1,figsize=(12,6))
ax.set_ylabel('Flux [normalised]')
ax.set_xlabel('Epoch [MJD]')
ax.set_title('Normalised ASASSN data {}'.format(fin))

ax.errorbar(tgout['MJD'],tgout['fnorm'],yerr=tgout['fnormerr'],label='g Band',fmt='.')
#ax.errorbar(tVout['MJD'],tVout['fnorm'],yerr=tVout['fnormerr'],label='V band',fmt='.')

ax.set_xlim([59000.,59500.])
ax.legend()

fig.savefig('_check_asassn9.pdf')

tn = vstack([tgout])
tn['Survey'] = "ASASSN"
tn['Source'] = obj
#plt.show()
tn.write(paths.data / 'obs_ASASSN-21co_ASASSN.ecsv',
    format='ascii.ecsv',overwrite=True)
