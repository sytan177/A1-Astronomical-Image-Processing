from astropy.io import fits
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

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

#(mu, sigma) = norm.fit(data)
mu, sigma = 3418.316, 11.364

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