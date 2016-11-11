import os
from astropy.io import fits
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
import numpy as np
import math
from math import sqrt, log


"""
BASIC REQUIREMENTS
"""

os.chdir('/Users/c1541417/Documents/891Paper/Maps')
arcsec = 3600   #Number of arcseconds in a degree
ngc891_coord = '35.639, 42.349'

"""
GET HEADER INFORMATION
"""

header_70 = fits.getheader('reprojected_image_70.fits')
#header_100 = fits.getheader('reprojected_image_100.fits')
#header_160 = fits.getheader('reprojected_image_160.fits')
#header_250 = fits.getheader('reprojected_image_250.fits')
#header_350 = fits.getheader('ngc891_spire_350.fits')

"""
READ IN THE DATA
"""

data_70 = fits.getdata('reprojected_image_70.fits')
#data_100 = fits.getdata('reprojected_image_100.fits')
#data_160 = fits.getdata('reprojected_image_160.fits')
#data_250 = fits.getdata('reprojected_image_250.fits')
#data_350 = fits.getdata('ngc891_spire_350.fits')
print 'works here'

"""
READ IN THE KERNEL FITS FILES
"""


kernel_70 = fits.getdata('reprojected_kernel.fits')
#kernel_100 = fits.getdata('Kernel_HiRes_PACS_100_to_SPIRE_350.fits')
#kernel_160 = fits.getdata('Kernel_HiRes_PACS_160_to_SPIRE_350.fits')
#kernel_250 = fits.getdata('Kernel_HiRes_SPIRE_250_to_SPIRE_350.fits')
print 'works here too'

"""
CONVOLVE
"""


convolved_70_to_350_new = convolve(data_70, kernel_70, boundary = 'wrap')
fits.writeto('70_convolved_to_350_test_new.fits', convolved_70_to_350_new, header_70, clobber=True)

#convolved_100_to_350_new = convolve(data_100, kernel_100, boundary='wrap')
#fits.writeto('100_convolved_to_350_test.fits', convolved_100_to_350_new, header_100, clobber=True)

#convolved_160_to_350_new = convolve(data_160, kernel_160, boundary='wrap')
#fits.writeto('160_convolved_to_350_test.fits', convolved_160_to_350_new, header_160, clobber=True)

#convolved_250_to_350_new = convolve(data_250, kernel_250, boundary='wrap')
#fits.writeto('250_convolved_to_350_test.fits', convolved_250_to_350_new, header_250, clobber=True)

