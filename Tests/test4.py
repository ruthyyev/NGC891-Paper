from astropy.io import fits 
import os
import numpy as np

os.chdir('/Users/c1541417/Documents/891Paper/Maps')

#data, header = fits.getdata('ngc891_spire_350_1.fits', header = True)

hdulist = fits.open('ngc891_spire_350_copy.fits')

hdu = hdulist[1]



new = np.delete(hdu.data, (0), axis = 0)

hdu.data = new

#fits.writeto('ngc891_spire_350_newnew.fits', hdu)

hdu.writeto('ngc891_spire_350_newnew.fits')