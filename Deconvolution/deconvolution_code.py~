from skimage import restoration
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import numpy as np



#where to find the maps:
os.chdir('/Users/c1541417/Documents/891Paper/Maps')


#import the PSF fits file:
psf = fits.getdata('spire_350_psf.fits')

psf = psf/psf.sum()

#fits.writeto('test_normalisation.fits', psf, clobber = False)

#import the image FITS file:
image, header = fits.getdata('ngc891_spire_350.fits', extname = 'image', header = True)
where_nan = np.where(np.isnan(image))
image[where_nan] = 0.


plt.imshow(psf)
plt.show()

plt.imshow(image)
plt.show()


#richardson_lucy(image, psf, iterations)
deconvolved = restoration.richardson_lucy(image, psf, 1)

plt.imshow(deconvolved)
plt.show()

fits.writeto('deconvolved_350_2.fits', deconvolved, clobber = True)
