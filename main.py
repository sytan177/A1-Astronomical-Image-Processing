# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from astropy.io import fits
# from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import numpy as np
# from scipy.ndimage.filters import generic_filter as gf
# from starcoord import star_coord
from numpy import unravel_index

# from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import chain

mu, sigma = 3418.316, 11.364

# masking_stars = star_coord
'''
def masking(centre, edge, coord):
    x = np.arange(0, 4611)
    y = np.arange(0, 2570)
    y0, x0 = int(centre[0]), int(centre[1])
    y1, x1 = int(edge[0]), int(edge[1])
    radius = int(np.sqrt((x1-x0)**2+(y1-y0)**2)) 
    mask = (y[np.newaxis,:]-y0)**2 + (x[:,np.newaxis]-x0)**2 < radius**2
    coord[mask] = 0.
    #coord[mask] = False

    return coord
'''


def find_source(data):
    m = max(data.ravel())
    if m > mu + 3*sigma:
        loc = np.argwhere(data.max() == data)
        return m, loc
    else:
        return 0

# data_selected = data_masked[1000:1100, 1000:1100]

def masking(centre, edge, coord):
    x = np.arange(0, 4611)
    y = np.arange(0, 2570)
    y0, x0 = int(centre[0]), int(centre[1])
    y1, x1 = int(edge[0]), int(edge[1])
    radius = int(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2))
    mask = (y[np.newaxis, :] - y0) ** 2 + (x[:, np.newaxis] - x0) ** 2 < radius ** 2
    coord[mask] = 0.
    # coord[mask] = False
    return coord


def find_rad(coord, loc):
    x0, y0 = loc[0]
    x, y = x0, y0
    while coord[x, y] > mu:
        x = x + 1
    radius = x - x0
    return radius


def local_bg_mask(centre, rad, annular_size=3):
    x0, y0 = centre[0]
    x = np.arange(0, 4611)
    y = np.arange(0, 2570)
    circle1 = (y[np.newaxis, :] - y0) ** 2 + (x[:, np.newaxis] - x0) ** 2 > rad ** 2

    circle2 = (y[np.newaxis, :] - y0) ** 2 + (x[:, np.newaxis] - x0) ** 2 < (rad + annular_size) ** 2

    circle = circle1 & circle2

    return circle


##### Background Masking ######
'''
hdulist = fits.open("A1_mosaic.fits")

image_data = hdulist[0].data
a = hdulist[0].data

hdu = hdulist[0]

hdulist.close()

coord = image_data

# data = [item for sublist in new_coord for item in sublist]

data_masked = np.zeros(image_data.shape)

data_masked[4415:4514, 114:2354] = coord[4415:4514, 114:2354]
data_masked[4212:4415, 114:2372] = coord[4212:4415, 114:2372]
data_masked[467:4212, 114:2471] = coord[467:4212, 114:2471]
data_masked[118:467, 116:1018] = coord[118:467, 116:1018]
data_masked[118:467, 1702:2471] = coord[118:467, 1702:2471]

data_masked[:, 1425:1448] = 0
data_masked[3203:3418, 770:780] = 0
data_masked[3707:3802, 2131:2137] = 0
data_masked[2704:2835, 968:976] = 0
data_masked[2224:2356, 901:908] = 0

while max(data_masked.ravel()) > 30000:
    m, centre = find_source(data_masked)
    print m, centre
    radius = find_rad(data_masked, centre)
    x0 = int(centre[0][0])
    y0 = int(centre[0][1])
    x = np.arange(0, 4611)
    y = np.arange(0, 2570)
    n = 2 * radius + 1
    mask_source = (y[np.newaxis, :] - y0) ** 2 + (x[:, np.newaxis] - x0) ** 2 < radius ** 2
    data_masked[mask_source] = 0.

hdulist[0].data = data_masked

try:
    hdulist.writeto('newtable.fits')
except IOError:
    os.remove('newtable.fits')
    hdulist.writeto('newtable.fits')

histogram = plt.hist(data_masked.ravel(), bins=20000)
plt.axis([3300, 3600, 0, 0.8e6])

plt.xlabel('Pixel Value')
plt.ylabel('Counts')
plt.show()

# best fit of data
hdulist = fits.open("newtable.fits")
image_data = hdulist[0].data
data = image_data.ravel()
data.sort()
i = 0
j = 0
data.sort()

while data[j] < 3350:
    j = j+1
while data[i] < 3480:
    i = i+1

data = data[j:i]
print len(data)
'''
#(mu, sigma) = norm.fit(data)
mu, sigma = 3418.316, 11.364
'''
# the histogram of the data
n, bins, patches = plt.hist(data, bins = 150, facecolor='green', alpha=0.75)
# add a 'best fit' line
y = 9076746*norm.pdf(bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)

#plot
plt.xlabel('Pixel value')
plt.ylabel('Counts')
plt.title(r'$\mathrm{Histogram\ of\ pixel:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)

plt.show()
'''
hdulist = fits.open("333photometryunder50002.fits")
data_masked = hdulist[0].data
'''
centre = [[1500,2000]]
fig = plt.figure()
ax = fig.add_subplot(111)
mask = local_bg(centre, 100)
ax.imshow(mask, aspect='auto', cmap=plt.cm.gray, interpolation='nearest')
'''
source_loc = []
source_rad = []
source_island = []    

while max(data_masked.ravel()) > mu + 3*sigma:
    m, centre = find_source(data_masked)
    print m
    radius = find_rad(data_masked, centre)
    if radius > 4:
        x0 = int(centre[0][0])
        y0 = int(centre[0][1])
        source_loc.append(centre)
        source_rad.append(radius)
        x = np.arange(0, 4611)
        y = np.arange(0, 2570)
        n = 2 * radius + 1

        mask_source = (y[np.newaxis,:]-y0)**2 + (x[:,np.newaxis]-x0)**2 < radius**2
        mask_bg = local_bg_mask(centre, radius)

        avg_bg_intensity = sum(data_masked[mask_bg]) / np.count_nonzero(mask_bg)
        raw_intensity = sum(data_masked[mask_source]) - avg_bg_intensity*np.count_nonzero(data_masked[mask_source])
        print raw_intensity
        source_island.append(raw_intensity)
        data_masked[mask_source] = 0.

    else:
        x0 = int(centre[0][0])
        y0 = int(centre[0][1])
        x = np.arange(0, 4611)
        y = np.arange(0, 2570)
        mask = (y[np.newaxis,:]-y0)**2 + (x[:,np.newaxis]-x0)**2 < 3**2
        data_masked[mask] = 0.

hdulist[0].data = data_masked

try:
    hdulist.writeto('333photometryunder34513sigma.fits' )
except IOError:
    os.remove('333photometryunder34513sigma.fits')
    hdulist.writeto('333photometryunder34513sigma.fits' )

with open('333source_location_34513sigma.txt', 'w') as f1:
    for item in source_loc:
        f1.write("%s\n" % item)
f1.close()

with open('333source_total_flux34513sigma.txt', 'w') as f2:
    for item in source_island:
        f2.write("%s\n" % item) 

f2.close()

with open('333source_total_rad34513sigma.txt', 'w') as f3:
    for item in source_rad:
        f3.write("%s\n" % item)      
f3.close()