import numpy as np
import scipy.ndimage
from astropy.io import fits
import os


#where to find the maps:
os.chdir('/Users/c1541417/Documents/891Paper/Maps')

data = fits.open('ngc891_spire_350.fits')
data2 = data[1]

psf = fits.open('spire_350_psf.fits')
psf2 = psf[1]


print 'Original array:'
print data2.data

factor =  data2.header['CDELT2'] / psf2.header['CDELT2']

print factor

zoomed = scipy.ndimage.zoom(data2.data, zoom = factor, order = 0)

print 'Resampled by "factor" with nearest interpolation:'
print zoomed

fits.writeto('ngc891_spire_350_new.fits', zoomed, data2.header, clobber = True)