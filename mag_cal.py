from astropy.io import fits
import numpy as np
import itertools
from methods import aper_size
from astropy.io import ascii
import matplotlib.pyplot as plt
import scipy.stats

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 6),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

loc_filename = '333source_location_3451.txt'
I_filename = '333source_total_flux3451.txt'

ZPinst = 2.530E+01 # Photometric ZP (mags) for default extinction

source_threshold = 3475.136
rad_threshold = aper_size()
dat_filename = "SDSS_source_counts_under_{}_aper_{}.dat".format(str(source_threshold), str(rad_threshold))

catalogue = ascii.read(dat_filename)
tol_counts = catalogue['TOL_COUNTS']

mag = []
err = []
mag_err = []

for i in range(len(tol_counts)):
    if tol_counts[i] <=  0:
        catalogue.remove_row(i)
    else:
        magi = -2.5*np.log10(float(tol_counts[i]))
        m = ZPinst + magi
        mag.append(m)
        err_in_counts = np.sqrt(tol_counts[i])
        err.append(err_in_counts)
        log_err = -2.5 * 0.434 * err_in_counts / tol_counts[i]
        mag_err.append(log_err)

catalogue['TOL_COUNTS_ERR'] = err
catalogue['MAG'] = mag
catalogue['MAG_ERR'] = mag_err

ascii.write(catalogue, dat_filename, overwrite = True)


### Plot magnitude and num of source fainter ###

def logNm(x):
    i = 0
    for j in range(len(mag)):
        if mag[j] <= x:
            i = i + 1

    return np.log10(i)

mag.sort()

y = []
for i in mag:
    y.append(logNm(i))

################# DATA OF THE FIRST SEGMENT ################
x1 = []
y1 = []
err_y1 = []

for j in range(len(tol_counts)):
    if mag[j] < 12.3:
        x1.append(mag[j])
        err_y1.append(mag_err[j])


for i in x1:
    y1.append(logNm(i))

k1, b1, r_value1, p_value1, std_err1 = scipy.stats.linregress(x1, y1)

fitting1 = []
for x_v in x1:
    fitting1.append(k1*x_v+b1)

################# DATA OF THE SECOND SEGMENT ################
x2 = []
y2 = []
err_y2 = []

for j in range(len(tol_counts)):
    if mag[j] > 12.3:
        if mag[j] < 17:
            x2.append(mag[j])
            err_y2.append(mag_err[j])

for i in x2:
    y2.append(logNm(i))

k2, b2, r_value2, p_value2, std_err2 = scipy.stats.linregress(x2, y2)

fitting2 = []
for x_v in x2:
    fitting2.append(k2*x_v+b2)

################# DATA OF THE THIRD SEGMENT ################
x3 = []
y3 = []
err_y3 = []

for j in range(len(tol_counts)):
    if mag[j] > 17:
        x3.append(mag[j])
        err_y3.append(mag_err[j])

for i in x3:
    y3.append(logNm(i))

k3, b3, r_value3, p_value3, std_err3 = scipy.stats.linregress(x3, y3)

fitting3 = []
for x_v in x3:
    fitting3.append(k3*x_v+b3)

'''
################### PLOT ######################
fig, axs = plt.subplots(nrows=4, ncols=1)
fig.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)

ax = axs[0]

ax.set_ylabel('$log_{10}N$' )
ax.errorbar(mag, y, yerr= mag_err , fmt='o', ms = 0.1)
ax.set_title('a)')

# With 4 subplots, reduce the number of axis ticks to avoid crowding.
ax.locator_params(nbins=4)

ax = axs[1]

ax.set_ylabel('$log_{10}N$')
ax.errorbar(x1, y1, yerr=err_y1, fmt='o', ms = 1)
ax.plot(x1, fitting1,label='y={:.2f}x{:.2f}'.format(k1,b1))
ax.text(0.5, 0.1, 'y={:.2f}x{:.2f}, $R^2$={:.2f}'.format(k1,b1,r_value1) ,size = 14, ha="center",
         transform=ax.transAxes)
ax.set_title('b) r < 12.3')

ax = axs[2]

ax.set_ylabel('$log_{10}N$')
ax.errorbar(x2, y2, yerr=err_y2, fmt='--o', ms = 1)
ax.plot(x2, fitting2)
ax.text(0.5, 0.1, 'y={:.2f}x{:.2f}, $R^2$={:.2f}'.format(k2,b2,r_value2) , size = 14,ha="center",
         transform=ax.transAxes)
ax.set_title('c) 12.3 < r < 17' )

ax = axs[3]
ax.set_xlabel('mag' )
ax.set_ylabel('$log_{10}N$' )
ax.errorbar(x3, y3, yerr=err_y3, fmt='--o', ms = 1)
ax.plot(x3, fitting3)
ax.text(0.5, 0.1, 'y={:.2f}x+{:.2f}, $R^2$={:.2f}'.format(k3,b3,r_value3) , size = 14,ha="center",
         transform=ax.transAxes)
ax.set_title('d) r > 17' )
'''
bin_num= 20
n, bins, patches = plt.hist(mag, bins = bin_num, facecolor='green', alpha=0.75)



plt.xlabel('mag')
plt.ylabel('$n_{gal}$')
plt.grid()


plt.show()