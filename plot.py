# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 09:43:45 2019

@author: St3916
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from methods import aper_size


source_threshold = 5000
rad_threshold = aper_size()
#dat_filename = "SDSS_source_counts_under_{}_aper_{}.dat".format(str(source_threshold), str(rad_threshold))
dat_filename = 'SDSS_source_counts_under_5000_aper_4.dat'

data = ascii.read(dat_filename)
mag =  data['MAG']

def num_of_source(x):
    i = 0
    for j in range(len(mag)):
        if mag[j] < x:
            i = i + 1
    return i
'''
x = np.linspace(10.5, 21, 20)
'''
y = []
for x_value in mag:
    Nm = num_of_source(x_value)
    y.append(np.log10(Nm))
    
k, b = np.polyfit(mag[:10], y[:10], 1)
print k, b

plt.scatter(mag, y)
plt.errorbar(mag, y, xerr=err, fmt='o')
plt.plot(mag, k*mag+b)
plt.show()

