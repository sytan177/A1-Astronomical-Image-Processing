from astropy.io import fits
import matplotlib.pyplot as plt
import os
import numpy as np

from methods import masking, find_source, find_rad, local_bg_mask

from astropy.io import ascii
from astropy.table import Table

mu, sigma = 3418.316, 11.364

hdulist = fits.open("333photometryunder50002.fits")
data_masked = hdulist[0].data

source_loc = []
source_rad = []
source_bg = []
source_island = []

source_threshold = 4999 # mu + 3*sigma:
rad_threshold = 4

while max(data_masked.ravel()) > source_threshold:
    m, centre = find_source(data_masked)
    print m
    radius = find_rad(data_masked, centre)
    
    if radius > rad_threshold:
        x0 = int(centre[0][0])
        y0 = int(centre[0][1])
        source_loc.append(centre)
        source_rad.append(radius)
        x = np.arange(0, 4611)
        y = np.arange(0, 2570)
        n = 2 * radius + 1

        mask_source = (y[np.newaxis,:]-y0)**2 + (x[:,np.newaxis]-x0)**2 < radius**2
        mask_bg = local_bg_mask(centre, radius)

        avg_bg_intensity = sum(data_masked[mask_bg]) / np.count_nonzero(data_masked(mask_bg))
        bg_to_remove = avg_bg_intensity*np.count_nonzero(data_masked[mask_source])
        source_bg.append(bg_to_remove)

        raw_intensity = sum(data_masked[mask_source]) - bg_to_remove
        print raw_intensity
        source_island.append(raw_intensity)
        data_masked[mask_source] = 0.

    else:
        x0 = int(centre[0][0])
        y0 = int(centre[0][1])
        x = np.arange(0, 4611)
        y = np.arange(0, 2570)
        mask = (y[np.newaxis,:]-y0)**2 + (x[:,np.newaxis]-x0)**2 < rad_threshold**2
        data_masked[mask] = 0.

hdulist[0].data = data_masked

fits_filename = "SDSS_source_counts_under_%s.fits" %source_threshold
dat_filename = "SDSS_source_counts_under_%s.dat" %source_threshold

try:
    hdulist.writeto(fits_filename)
except IOError:
    os.remove(fits_filename)
    hdulist.writeto(fits_filename)

source_to_store = Table([source_loc, source_island, source_bg, source_rad], names = ['XY', 'TOL_COUNTS', 'TOL_BG', 'APER_RAD'])

ascii.write(source_to_store, dat_filename)