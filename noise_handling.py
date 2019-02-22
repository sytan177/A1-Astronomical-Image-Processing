from astropy.io import fits
import cv2
#from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf
#from starcoord import star_coord
from numpy import unravel_index

#from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import chain

from astropy.visualization import HistEqStretch
from astropy.visualization.mpl_normalize import ImageNormalize

mu_i, sigma_i = 3418.316, 11.364

hdulist = fits.open("333photometryunder5000.fits")
image_data = hdulist[0].data
hdulist.close()

image_data = gf(image_data, sigma = sigma_i) - mu_i


norm = ImageNormalize(stretch=HistEqStretch(image_data))

hdulist[0].data = image_data



try:
    hdulist.writeto('newtable1.fits' )
except IOError:
    os.remove('newtable1.fits')
    hdulist.writeto('newtable1.fits' )