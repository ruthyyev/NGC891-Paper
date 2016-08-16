"""
A code to convolve images to the lowest resolution image after images 
have been regridded using ~/Documents/891Paper/reprojection_code.py.
"""

import os
from astropy.io import fits
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
import numpy as np
import math
from math import sqrt, log

#------------------------------------------------------------------------------
"""BASIC REQUIREMENTS"""
#------------------------------------------------------------------------------

os.chdir('/Users/c1541417/Documents/891Paper/Maps')
arcsec = 3600   #Number of arcseconds in a degree
ngc891_coord = '35.639, 42.349'

#------------------------------------------------------------------------------
"""READ IN FITS FILES"""
#------------------------------------------------------------------------------

#Get header information:
header_1mm = fits.getheader('reprojected_image_1mm.fits')
header_2mm = fits.getheader('reprojected_image_2mm.fits')
header_70 = fits.getheader('reprojected_image_70.fits')
header_100 = fits.getheader('reprojected_image_100.fits')
header_160 = fits.getheader('reprojected_image_160.fits')
header_250 = fits.getheader('reprojected_image_250.fits')
header_350 = fits.getheader('ngc891_spire_350.fits')

#Get the pixel size (in arcsec) from the header files:
pix_size_1mm = header_1mm['CDELT2']*arcsec
pix_size_2mm = header_2mm['CDELT2']*arcsec
pix_size_70  = header_70['CDELT2']*arcsec
pix_size_100 = header_100['CDELT2']*arcsec
pix_size_160 = header_160['CDELT2']*arcsec
pix_size_250 = header_250['CDELT2']*arcsec
pix_size_350 = header_350['CDELT2']*arcsec

#Make each pixel size from each band into a numpy array:
pix_size = np.array([pix_size_1mm, pix_size_2mm, pix_size_70, pix_size_100, pix_size_160, pix_size_250, pix_size_350])

#------------------------------------------------------------------------------
"""DEFINE THE POINT SPREAD FUNCTIONS/BEAM WIDTHS"""
#------------------------------------------------------------------------------                

"""

Resulting values have units of pixels. 

The average PSF for PACS observations were taken from the 
PACS observers manual, Table 3.1 which can be found at:
http://herschel.esac.esa.int/Docs/PACS/pdf/pacs_om.pdf
We used values at 60 arcsec/s scanning speed.

The average beam widths for SPIRE observations were taken from:
http://www.aanda.org/articles/aa/pdf/2010/10/aa14519-10.pdf

NIKA2 beam widths are taken from: 
http://www.iram.fr/wiki/nika2/images/4/49/Commissioning.pdf

"""

#NIKA 1 and 2 mm bands:
av_beam_1mm_pix = 12/pix_size[0]

av_beam_2mm_pix = 18/pix_size[1]


#Herschel PACS bands:
av_PSF_70 = (5.75+9)/2
av_PSF_70_pix = av_PSF_70/pix_size[2]      

av_PSF_100 = (6.89+9.74)/2
av_PSF_100_pix = av_PSF_100/pix_size[3]

av_PSF_160 = (11.31+13.32)/2
av_PSF_160_pix =av_PSF_160/pix_size[4] 
    

#Herschel SPIRE bands:
av_beam_250_pix = 18.1/pix_size[5]  
            
av_beam_350_pix = 24.9/pix_size[6]             


#The amount needed to convolve each beam to match the 350 micron beam:
beam_1mm_to_350 = sqrt(((av_beam_350_pix)**2)-((av_beam_1mm_pix)**2))
beam_2mm_to_350 = sqrt((av_beam_350_pix**2)-(av_beam_2mm_pix**2))
beam_70_to_350 = sqrt((av_beam_350_pix**2)-(av_PSF_70_pix**2))
beam_100_to_350 = sqrt((av_beam_350_pix**2)-(av_PSF_100_pix**2))
beam_160_to_350 = sqrt((av_beam_350_pix**2)-(av_PSF_160_pix**2))
beam_250_to_350 = sqrt((av_beam_350_pix**2)-(av_beam_250_pix**2))


##These four lines of code print the values of the beam change in 
#arcseconds, i.e. how much should be added to the 9.01 arcsecs beam to 
#convolve it to the 13.645arcsec beam. Same for the other three:
print 'beam_1mm_to_350 = ' +str(beam_1mm_to_350)
print 'beam_2mm_to_350 = ' +str(beam_2mm_to_350)
print 'beam_70_to_350 = '+ str(beam_70_to_350)
print 'beam_100_to_350 = ' + str(beam_100_to_350)                                                                                 
print 'beam_160_to_350 = ' + str(beam_160_to_350)
print 'beam_250_to_350  = ' + str(beam_250_to_350)

#------------------------------------------------------------------------------
"""CALCULATE SIGMA"""
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
sigma_1mm_to_350 = beam_1mm_to_350/const
sigma_2mm_to_350 = beam_2mm_to_350/const
sigma_70_to_350 = beam_70_to_350/const      
sigma_100_to_350 = beam_100_to_350/const 
sigma_160_to_350 = beam_160_to_350/const
sigma_250_to_350 = beam_250_to_350/const


#Print the values of sigma:
print 'sigma_1mm_to_350 = ' +str(sigma_1mm_to_350)
print 'sigma_2mm_to_350 = ' +str(sigma_2mm_to_350)
print 'sigma_70_to_350 = '+ str(sigma_70_to_350)
print 'sigma_100_to_350 = ' + str(sigma_100_to_350)                                                                                 
print 'sigma_160_to_350 = ' + str(sigma_160_to_350)
print 'sigma_250_to_350  = ' + str(sigma_250_to_350)





"""

Here follow the conversion factors for the 70um, 160um, 250um, 350um and 500um 
from MJy/sr to Jy/pixel. It has been calculated as follows: 

"""

#------------------------------------------------------------------------------
"""CONVERT MJY/SR TO JY/PIXEL"""
#------------------------------------------------------------------------------


"""

1 radian = 180/pi degrees 
1 steradian = (180/pi) square degrees
1 deg = 3600 arcsecs hence 
180 deg = 180*3600 arcsecs = 648000 arcsecs
So 1 steradian = (648000^2/3.2^2	) square arcsecs = 41006249999.99999 square arcsecs or 
1 square arcsec = 2.4386526444139618e-11 steradian
So Mjy/sr = 10^6Jy/sr * 1sr/0.4253*10^11 square arcsec = 10^6Jy/41006249999.99999 square arcsec = 2.4386526444139618*10^(-5) Jy/square arcsec

1 sr = 41006249999.99999 arcsec^2
1 arcsec^2 = 2.4386526444139618e-11 sr
1 MJy/sr = 2.4386526444139618*10^(-5) Jy/arcsec^2 ---> this is called 'x' below

"""

x = 2.4386526444139618e-5



#NIKA2 1mm band:
pix1 = pix_size_1mm**2
factor_1mm = 1/pix1  
jy_per_pix_1mm = factor_1mm/x


#NIKA2 2mm band:
pix2 = pix_size_2mm**2
factor_2mm = 1/pix2 
jy_per_pix_2mm = factor_2mm/x


#PACS 70 band:
pix3 = pix_size_70**2
factor_70 = 1/pix3
jy_per_pix_70 = factor_70/x


#PACS 100 band:
pix4 = pix_size_100**2
factor_100 = 1/pix4 
jy_per_pix_100 = factor_100/x


#PACS 160 band:
pix5 = pix_size_160**2
factor_160 = 1/pix5 
jy_per_pix_160 = factor_160/x


#SPIRE 250 band:
pix6 = pix_size_250**2
factor_250 = 1/pix6 
jy_per_pix_250 = factor_250/x


#SPIRE 350 band:
pix7 = pix_size_350**2
factor_350 = 1/pix7 
jy_per_pix_350 = factor_350/x

"""

Get the numpy data for the original images. We will convert this to 
Jy/pixel and then the conversion should be carried over to the rest of the 
processes.

"""

data_1mm = fits.getdata('reprojected_image_1mm.fits')

data_2mm = fits.getdata('reprojected_image_2mm.fits')

data_70 = fits.getdata('reprojected_image_70.fits')

data_100 = fits.getdata('reprojected_image_100.fits')

data_160 = fits.getdata('reprojected_image_160.fits')

data_250 = fits.getdata('reprojected_image_250.fits')

data_350 = fits.getdata('ngc891_spire_350.fits')


#Correct the data by multiplying with the conversion factor for the 1mm:
data_1mm_corrected = (data_1mm)*factor_1mm
fits.writeto('ngc891_nika2_1mm_corrected.fits', data_1mm_corrected, header_1mm, clobber=True)

data_2mm_corrected = (data_2mm)*factor_2mm
fits.writeto('ngc891_nika2_2mm_corrected.fits', data_2mm_corrected, header_2mm, clobber=True)

data_70_corrected = (data_70)*factor_70
fits.writeto('ngc891_nika2_70_corrected.fits', data_1mm_corrected, header_70, clobber=True)

data_100_corrected = (data_100)*factor_100
fits.writeto('ngc891_nika2_100_corrected.fits', data_100_corrected, header_100, clobber=True)

data_160_corrected = (data_160)*factor_160
fits.writeto('ngc891_nika2_160_corrected.fits', data_160_corrected, header_160, clobber=True)

data_250_corrected = (data_250)*factor_250
fits.writeto('ngc891_nika2_250_corrected.fits', data_160_corrected, header_250, clobber=True)

#------------------------------------------------------------------------------
"""CONVOLVE"""
#------------------------------------------------------------------------------

"""

Here, the Gaussian2DKernel is defined (imported from astropy.convolution) and 
used. The parameter used is the standard deviation, sigma, in pixels for each 
of the different wavelengths. 

"""

gauss_1mm_to_350 = Gaussian2DKernel(sigma_1mm_to_350)
gauss_2mm_to_350 = Gaussian2DKernel(sigma_2mm_to_350)
gauss_70_to_350 = Gaussian2DKernel(sigma_70_to_350)  
gauss_100_to_350 = Gaussian2DKernel(sigma_100_to_350)  
gauss_160_to_350 = Gaussian2DKernel(sigma_160_to_350)
gauss_250_to_350 = Gaussian2DKernel(sigma_250_to_350)

#Convolve:
convolved_1mm_to_350_new = convolve(data_1mm_corrected, gauss_1mm_to_350, boundary='wrap')
fits.writeto('1mm_convolved_to_350.fits', convolved_1mm_to_350_new, header_1mm, clobber=True)

convolved_2mm_to_350_new = convolve(data_2mm_corrected, gauss_2mm_to_350, boundary='wrap')
fits.writeto('2mm_convolved_to_350.fits', convolved_2mm_to_350_new, header_2mm, clobber=True)

convolved_70_to_350_new = convolve(data_70_corrected, gauss_70_to_350, boundary='wrap')
fits.writeto('70_convolved_to_350.fits', convolved_70_to_350_new, header_70, clobber=True)

convolved_100_to_350_new = convolve(data_100_corrected, gauss_100_to_350, boundary='wrap')
fits.writeto('100_convolved_to_350.fits', convolved_100_to_350_new, header_100, clobber=True)

convolved_160_to_350_new = convolve(data_160_corrected, gauss_160_to_350, boundary='wrap')
fits.writeto('160_convolved_to_350.fits', convolved_160_to_350_new, header_160, clobber=True)

convolved_250_to_350_new = convolve(data_250_corrected, gauss_250_to_350, boundary='wrap')
fits.writeto('250_convolved_to_350.fits', convolved_250_to_350_new, header_250, clobber=True)






