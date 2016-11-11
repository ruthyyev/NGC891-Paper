from astropy.io import fits
import os

os.chdir('/Users/c1541417/Documents/891Paper/Maps')

data, header = fits.getdata("Kernel_HiRes_PACS_70_to_SPIRE_350.fits", header=True)

header['CRVAL1']  =    35.64035037478851  
header['CRVAL2']  =    42.35599876715473
header['CDELT2']  =    0.0000694444444444444

fits.writeto('Kernel_HiRes_PACS_70_to_SPIRE_350_new.fits', data, header, clobber=True)