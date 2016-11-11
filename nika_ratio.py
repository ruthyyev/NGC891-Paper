"""A code to regrid NIKA2 1mm map to the 2mm, convolve with a gaussian 
and find the ratio of the images"""

#------------------------------------------------------------------------------
#CALL MODULES
#------------------------------------------------------------------------------

from astropy.io import fits
import numpy as np
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
#from astropy.utils.data import get_pkg_data_filename
#from astropy.wcs import WCS
import matplotlib.pyplot as plt
import os
from reproject import reproject_exact
import math
from math import sqrt, log

#------------------------------------------------------------------------------
#BASIC REQUIREMENTS
#------------------------------------------------------------------------------

os.chdir('/Users/c1541417/Documents/891Paper/Maps')
arcsec = 3600   #Number of arcseconds in a degree
ngc891_coord = '35.639, 42.349'

#------------------------------------------------------------------------------
#IMPORT FITS IMAGES
#------------------------------------------------------------------------------

hdulist1 = fits.open('ngc891_nika2_1mm.fits')
hdu1mm = hdulist1[1]

hdulist2 = fits.open('ngc891_nika2_2mm.fits')
hdu2mm = hdulist2[1]


#------------------------------------------------------------------------------
#VIEW IMAGES TO TEST FITS HANDLING
#------------------------------------------------------------------------------

"""
plt.imshow(hdu1mm.data, origin = 'lower')
plt.show()

plt.imshow(hdu2mm.data, origin = 'lower')
plt.show()
"""


#------------------------------------------------------------------------------
#REPROJECT IMAGES
#------------------------------------------------------------------------------

new_image_1mm, footprint = reproject_exact(hdu1mm, hdu2mm.header)

#------------------------------------------------------------------------------
#SAVE REPROJECTED IMAGES
#------------------------------------------------------------------------------

fits.writeto('reprojected_image_1mm_to_2mm.fits', new_image_1mm, hdu2mm.header, clobber = True)

#------------------------------------------------------------------------------
#FIND THE BEAM RATIOS TO CONVOLVE WITH
#------------------------------------------------------------------------------

#Import header files
header_1mm = fits.getheader('reprojected_image_1mm_to_2mm.fits')
header_2mm = fits.getheader('ngc891_nika2_2mm.fits', 1)

#Get the pixel size (in arcsec) from the header files:
pix_size_1mm = header_1mm['CDELT2']*arcsec
pix_size_2mm = header_2mm['CDELT2']*arcsec

#Make each pixel size from each band into a numpy array:
pix_size = np.array([pix_size_1mm, pix_size_2mm])

#Define the PSFs:
av_beam_1mm_pix = 12/pix_size[0]

av_beam_2mm_pix = 18/pix_size[1]

#The amount needed to convolve each beam to match the 2mm beam:
beam_1mm_to_2mm = sqrt(((av_beam_2mm_pix)**2)-((av_beam_1mm_pix)**2))

#Print the values of the beam change in arcseconds, i.e. how much should be 
#added to the beam to convolve it to the larger beam
print 'beam_1mm_to_2mm = ' +str(beam_1mm_to_2mm)


#------------------------------------------------------------------------------
#CALCULATE SIGMA
#------------------------------------------------------------------------------

"""
The equation for sigma, the standard deviation in arcseconds, given the beam 
at FWHM is : sigma = FWHM/2*sqrt(2*ln2). In this step the denominator is named 
const to simplify the code.
"""

#Standard deviation in arcsec:
const=2*sqrt(2*log(2))    
print 'const = ',float(const)                                


#Express sigma in pixel size as opposed to arcsec:
sigma_1mm_to_2mm = beam_1mm_to_2mm/const

#Print the value of sigma:
print 'sigma_1mm_to_2mm = ' +str(sigma_1mm_to_2mm)

gauss_1mm_to_2mm = Gaussian2DKernel(sigma_1mm_to_2mm)


#------------------------------------------------------------------------------
#CONVOLVE
#------------------------------------------------------------------------------

convolved_1mm_to_2mm_new = convolve(new_image_1mm, gauss_1mm_to_2mm, boundary='wrap')
fits.writeto('1mm_convolved_to_2mm_new.fits', convolved_1mm_to_2mm_new, header_1mm, clobber=True)

#------------------------------------------------------------------------------
#NORMALISE AND FIND THE RATIOS
#------------------------------------------------------------------------------

norm_1mm = convolved_1mm_to_2mm_new / convolved_1mm_to_2mm_new.sum()
norm_2mm = hdu2mm.data / hdu2mm.data.sum()

fits.writeto('test_normalisation2.fits', norm_2mm, clobber = False)

#------------------------------------------------------------------------------
#PLOT AND SAVE THE RATIO MAPS
#------------------------------------------------------------------------------
"""
ratio_map1 = hdu2mm.data / convolved_1mm_to_2mm_new
fits.writeto('ratio_map.fits', ratio_map1, clobber = True)

ratio_map2 = convolved_1mm_to_2mm_new / hdu2mm.data
fits.writeto('ratio_map2.fits', ratio_map2, clobber = True)
"""

ratio_map1 = norm_2mm / norm_1mm
fits.writeto('ratio_map_norm.fits', ratio_map1, clobber = True)

ratio_map2 = norm_1mm / norm_2mm
fits.writeto('ratio_map2_norm.fits', ratio_map2, clobber = True)



 


