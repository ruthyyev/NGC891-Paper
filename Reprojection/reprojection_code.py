"""A code to regrid NIKA2 and Herschel pixels to 250 micron SPIRE pixels"""

#------------------------------------------------------------------------------
"""CALL MODULES"""
#------------------------------------------------------------------------------

from astropy.io import fits
import numpy as np
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

hdulist2 = fits.open('ngc891_nika2_1mm.fits')
hdu1mm = hdulist2[1]

hdulist3 = fits.open('ngc891_nika2_2mm.fits')
hdu2mm = hdulist3[1]

hdulist4 = fits.open('ngc891_pacs_70.fits')
hdu70 = hdulist4[0]

hdulist5 = fits.open('ngc891_pacs_100.fits')
hdu100 = hdulist5[0]

hdulist6 = fits.open('ngc891_pacs_160.fits')
hdu160 = hdulist6[0]

hdulist7 = fits.open('ngc891_spire_250.fits')
hdu250 = hdulist7[1]


#------------------------------------------------------------------------------
"""VIEW IMAGES TO TEST FITS HANDLING"""
#------------------------------------------------------------------------------

plt.imshow(hdu350.data, origin = 'lower')
plt.show()

plt.imshow(hdu2mm.data, origin = 'lower')
plt.show()

plt.imshow(hdu70.data, origin = 'lower')
plt.show()



#------------------------------------------------------------------------------
"""REPROJECT IMAGES"""
#------------------------------------------------------------------------------

new_image_1mm, footprint = reproject_exact(hdu1mm, hdu350.header)

new_image_2mm, footprint = reproject_exact(hdu2mm, hdu350.header)

new_image_70, footprint = reproject_exact(hdu70, hdu350.header)

new_image_100, footprint = reproject_exact(hdu100, hdu350.header)

new_image_160, footprint = reproject_exact(hdu160, hdu350.header)

new_image_250, footprint = reproject_exact(hdu250, hdu350.header)


#------------------------------------------------------------------------------
"""SAVE REPROJECTED IMAGES"""
#------------------------------------------------------------------------------

fits.writeto('reprojected_image_1mm.fits', new_image_1mm, hdu350.header)

fits.writeto('reprojected_image_2mm.fits', new_image_2mm, hdu350.header) 

fits.writeto('reprojected_image_70.fits', new_image_70, hdu350.header, clobber = True) 

fits.writeto('reprojected_image_100.fits', new_image_100, hdu350.header) 

fits.writeto('reprojected_image_160.fits', new_image_160, hdu350.header) 

fits.writeto('reprojected_image_250.fits', new_image_250, hdu350.header) 



 


