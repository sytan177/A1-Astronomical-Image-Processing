from astropy.io import fits
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
from methods import masking, find_source, find_rad, local_bg_mask

############ Background Masking #############

hdulist = fits.open("A1_mosaic.fits")
image_data = hdulist[0].data
coord = image_data
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
    #print m, centre
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
