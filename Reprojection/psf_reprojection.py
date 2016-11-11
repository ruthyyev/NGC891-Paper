"""A code to match a map with its PSF"""

#------------------------------------------------------------------------------
"""CALL MODULES"""
#------------------------------------------------------------------------------

from astropy.io import fits
#import numpy as np
#from astropy.utils.data import get_pkg_data_filename
#from astropy.wcs import WCS
import matplotlib.pyplot as plt
import os
from reproject import reproject_exact

#------------------------------------------------------------------------------
"""BASIC REQUIREMENTS"""
#------------------------------------------------------------------------------

os.chdir('/Users/c1541417/Documents/891Paper/Maps')

#------------------------------------------------------------------------------
"""IMPORT FITS IMAGES"""
#------------------------------------------------------------------------------

hdulist1 = fits.open('ngc891_spire_350.fits')
hdu350 = hdulist1[1]

hdulist2 = fits.open('spire_350_psf_new.fits')
hdupsf = hdulist2[0]



plt.imshow(hdu350.data)
plt.show()
plt.imshow(hdupsf.data)
plt.show()



#------------------------------------------------------------------------------
"""REPROJECT IMAGES"""
#------------------------------------------------------------------------------

#new_image_350, footprint = reproject_exact(hdu350, hdupsf.header, shape_out = (229, 229), hdu_in = hdu350)
new_spire_350, footprint = reproject_exact(hdu350, hdupsf.header)


plt.imshow(new_spire_350)
plt.show()

#------------------------------------------------------------------------------
"""SAVE REPROJECTED IMAGES"""
#------------------------------------------------------------------------------

#fits.writeto('new_spire_350.fits', new_image_350, hdupsf.header, clobber = True)
fits.writeto('snew_spirepire_350.fits', new_spire_350, hdupsf.header, clobber = True) 

#plt.plot(footprint)
#plt.show()



"""
//anaconda/lib/python2.7/site-packages/numpy/ma/core.py:4085: UserWarning: Warning: converting a masked element to nan.
  warnings.warn("Warning: converting a masked element to nan.")
"""

 


