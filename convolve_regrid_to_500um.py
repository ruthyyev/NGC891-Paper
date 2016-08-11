#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""This is a code to a) regrid an image from one pixel size to the pixel size of 
another, target image, b) convolve the regridded image to a specified, new beam 
size determined  by the beam sizes of the #regridded image and the target image 
and c) perfom a division between the regridded image and the target image. The 
aim of this division is to obtain ratio intensity maps in an attempt to observe 
any dust in the region of interest. The fits files used here are Hi-GAL files. 
The Hi-GAL data was obtained in parallel scan mode. (Traficante et al, 2011, 
"Data reduction pipeline for the Hi-GAL survey").

In the convolution part of the code, the argument in the Gaussian2DKernel 
(see astropy documentation, 
http://docs.astropy.org/en/stable/api/astropy.convolution.Gaussian2DKernel.html),  
is the standard #deviation in pixels at FWHM. To get this value, given the 
beam width, the approach would be to find the standard deviation (sigma) 
in arcseconds and then, for each wavelength band, divide that value with 
the pixel size to get sigma in terms of pixels instead of arcseconds."""


#------------------------------------------------------------------------------
"""REQUIRED PACKAGES"""              
#------------------------------------------------------------------------------

import astropy
import os
import aplpy
from astropy.io import fits
import numpy
import montage_wrapper
from montage_wrapper import mHdr,mProject
from astropy.convolution import convolve, Gaussian2DKernel
import numpy as np
from numpy import array
import math
from math import sqrt, log                            
import cmath
#import montage_wrapper
#from montage_wrapper import mProject, mHdr


#------------------------------------------------------------------------------
"""BASIC REQUIREMENTS"""
#------------------------------------------------------------------------------                    


os.chdir('/Users/c1541417/Documents/891Paper/Maps')
arcsec = 3600   #Number of arcseconds in a degree
ngc891_coord = '35.639, 42.349'    
print 'ngc891_coord = ' + ngc891_coord           
width = 0.5     #Size of the grid in degrees.


#------------------------------------------------------------------------------
"""READ IN FITS FILES"""
#------------------------------------------------------------------------------


"""Read in the headers from the HIGAL fits files to be used here and extract the 
CDELT values in order to calculate the pixel sizes in arcseconds."""



#Get header information:
header_1mm = fits.getheader('ngc891_nika2_1mm.fits', 1)
header_2mm = fits.getheader('ngc891_nika2_2mm.fits', 1)
header_70 = fits.getheader('ngc891_pacs_70.fits')
header_100 = fits.getheader('ngc891_pacs_100.fits')
header_160 = fits.getheader('ngc891_pacs_160.fits')
header_250 = fits.getheader('ngc891_spire_250.fits')
header_350 = fits.getheader('ngc891_spire_350.fits')


#Get the pixel size (in arcsec) from the header files:
pix_size_1mm = header_1mm['CDELT2']*arcsec
pix_size_2mm = header_2mm['CDELT2']*arcsec
pix_size_70  = header_70['CDELT2']*arcsec
pix_size_100 = header_100['CDELT2']*arcsec
pix_size_160 = header_160['CDELT2']*arcsec
pix_size_250 = header_250['CDELT2']*arcsec
pix_size_350 = header_350['CDELT2']*arcsec


print pix_size_1mm
print pix_size_2mm
print pix_size_70
print pix_size_100
print pix_size_160
print pix_size_250
print pix_size_350


#Make each pixel size from each band into a numpy array:
pix_size = np.array([pix_size_1mm, pix_size_2mm, pix_size_70, pix_size_100, pix_size_160, pix_size_250, pix_size_350])


#------------------------------------------------------------------------------
"""PIXEL SIZE RATIOS"""
#------------------------------------------------------------------------------                


"""Here follow the area ratios between the different pixel areas. These will be 
used later in order to correct for the effect the regridding has on surface 
brightness."""

pix_ratio_350_by_1mm = (pix_size[6]/pix_size[0]**2)
pix_ratio_350_by_2mm = (pix_size[6]/pix_size[1]**2)
pix_ratio_350_by_70 = (pix_size[6]/pix_size[2]**2)
pix_ratio_350_by_100 = (pix_size[6]/pix_size[3]**2)
pix_ratio_350_by_160 = (pix_size[6]/pix_size[4]**2)
pix_ratio_350_by_250 = (pix_size[6]/pix_size[5]**2)


print 'pix_ratio_350_by_1mm = ' + str(pix_ratio_350_by_1mm)
print 'pix_ratio_350_by_2mm = ' + str(pix_ratio_350_by_2mm)
print 'pix_ratio_350_by_70 = ' + str(pix_ratio_350_by_70)
print 'pix_ratio_350_by_100 = ' + str(pix_ratio_350_by_100)
print 'pix_ratio_350_by_160 = ' + str(pix_ratio_350_by_160)
print 'pix_ratio_350_by_250 = ' + str(pix_ratio_350_by_250)


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



#Print the number of pixels per beam:
print 'av_beam_1mm_pix = ' +str(av_beam_1mm_pix)
print 'av_beam_2mm_pix = ' +str(av_beam_2mm_pix)
print 'av_PSF_70_pix = ' + str(av_PSF_70_pix)
print 'av_PSF_100_pix = ' +str(av_PSF_100_pix)
print 'av_PSF_160_pix = ' + str(av_PSF_160_pix)
print 'av_beam_250_pix = ' + str(av_beam_250_pix)
print 'av_beam_350_pix = ' + str(av_beam_350_pix)


#The amount needed to convolve each beam to match the 350 micron beam:
beam_1mm_to_350 = sqrt(((av_beam_350_pix)**2)-((av_beam_1mm_pix)**2))
beam_2mm_to_350 = sqrt((av_beam_350_pix**2)-(av_beam_2mm_pix**2))
beam_70_to_350 = sqrt((av_PSF_70_pix**2)-(av_beam_350_pix**2))
beam_100_to_350 = sqrt((av_PSF_100_pix**2)-(av_beam_350_pix**2))
beam_160_to_350 = sqrt((av_PSF_160_pix**2)-(av_beam_350_pix**2))
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


#------------------------------------------------------------------------------
"""CREATE HEADER WITH CO-ORDINATES OF THE GALAXY"""
#------------------------------------------------------------------------------


"""

This section of the code creates a header with central coordinates those of 
the galaxy, grid_size big enough to include the remnant and pixel size which 
depends upon which wavelength is chosen each time. Seven headers are created, 
each for each wavelength. 

"""


target_header_1mm = mHdr(ngc891_coord, width, 'target_1mm.txt', pix_size=pix_size[0])
target_header_2mm = mHdr(ngc891_coord, width, 'target_2mm.txt', pix_size=pix_size[1])
target_header_70 = mHdr(ngc891_coord, width, 'target_70.txt', pix_size=pix_size[2])
target_header_100 = mHdr(ngc891_coord, width, 'target_100.txt', pix_size=pix_size[3])
target_header_160 = mHdr(ngc891_coord, width, 'target_160.txt', pix_size=pix_size[4])            
target_header_250 = mHdr(ngc891_coord, width, 'target_250.txt', pix_size=pix_size[5])
target_header_350 = mHdr(ngc891_coord, width, 'target_350.txt', pix_size=pix_size[6])



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



#------------------------------------------------------------------------------
"""CONVERT AND REGRID THE ORIGINAL IMAGES"""
#------------------------------------------------------------------------------

"""

Get the numpy data for the original images. We will convert this to 
Jy/pixel and then the conversion should be carried over to the rest of the 
processes.

"""


data_1mm = fits.getdata('ngc891_nika2_1mm.fits')
print(data_1mm)


data_2mm = fits.getdata('ngc891_nika2_2mm.fits')
print(data_2mm)

data_70 = fits.getdata('ngc891_pacs_70.fits')
print(data_70)

data_100 = fits.getdata('ngc891_pacs_100.fits')
print(data_100)

data_160 = fits.getdata('ngc891_pacs_160.fits')
print(data_160)

data_250 = fits.getdata('ngc891_spire_250.fits')
print(data_250)

data_350 = fits.getdata('ngc891_spire_350.fits')
print(data_350)


#Get the header to use in order to save the corrected file in fits format:
header_1mm = fits.getheader('ngc891_nika2_1mm.fits', 1)
header_2mm = fits.getheader('ngc891_nika2_2mm.fits', 1)
header_70 = fits.getheader('ngc891_pacs_70.fits')
header_100 = fits.getheader('ngc891_pacs_100.fits')
header_160 = fits.getheader('ngc891_pacs_160.fits')
header_250 = fits.getheader('ngc891_spire_250.fits')
header_350 = fits.getheader('ngc891_spire_350.fits')

#------------------------------------------------------------------------------
"""NIKA2 1MM IMAGES - CONVOLVE/REGRID"""
#------------------------------------------------------------------------------


#Correct the data by multiplying with the conversion factor for the 1mm:
data_1mm_corrected = (data_1mm)*factor_1mm
print(data_1mm_corrected)
fits.writeto('ngc891_nika2_1mm_corrected.fits', data_1mm_corrected, header_1mm)


#Convolve:
convolved_image_1mm_to_350_new = convolve(data_1mm_corrected, gauss_1mm_to_350)
fits.writeto('1mm_convolved_to_350.fits', convolved_image_1mm_to_350_new, header_1mm)


#Regrid:
regridded_1mm = mProject('ngc891_nika2_1mm_corrected.fits', '1mm_regridded.fits','target_1mm.txt')

regridded_1mm_to_350 = mProject('1mm_convolved_to_350.fits', '1mm_regridded_to_350.fits', 'target_350.txt')

data_regridded_1mm_to_350 = fits.getdata('1mm_regridded_to_350.fits')

header_regridded_1mm_to_350 = fits.getheader('1mm_regridded_to_350.fits')

data_regridded_1mm_to_350_scaled = data_regridded_1mm_to_350*pix_ratio_350_by_1mm

fits.writeto('1mm_regridded_to_350.fits', data_regridded_1mm_to_350_scaled, header_regridded_1mm_to_350, clobber=True)



#------------------------------------------------------------------------------
"""NIKA2 2MM IMAGES - CONVOLVE/REGRID"""
#------------------------------------------------------------------------------


#Correct the data by multiplying with the conversion factor for the 2mm:
data_2mm_corrected = (data_2mm)*factor_2mm
print(data_2mm_corrected)
fits.writeto('ngc891_nika2_2mm_corrected.fits', data_2mm_corrected, header_2mm)


#Convolve:
convolved_image_2mm_to_350_new = convolve(data_2mm_corrected, gauss_2mm_to_350)
fits.writeto('2mm_convolved_to_350.fits', convolved_image_2mm_to_350_new, header_2mm)


#Regrid:
regridded_2mm = mProject('ngc891_nika2_2mm_corrected.fits', '2mm_regridded.fits','target_2mm.txt')

regridded_2mm_to_350 = mProject('2mm_convolved_to_350.fits', '2mm_regridded_to_350.fits', 'target_350.txt')

data_regridded_2mm_to_350 = fits.getdata('2mm_regridded_to_350.fits')

header_regridded_2mm_to_350 = fits.getheader('2mm_regridded_to_350.fits')

data_regridded_2mm_to_350_scaled = data_regridded_2mm_to_350*pix_ratio_350_by_2mm

fits.writeto('2mm_regridded_to_350.fits', data_regridded_2mm_to_350_scaled, header_regridded_2mm_to_350, clobber=True)



#------------------------------------------------------------------------------
"""PACS 70 MICRON IMAGES - CONVOLVE/REGRID"""
#------------------------------------------------------------------------------


#Correct the data by multiplying with the conversion factor for the 70:
data_70_corrected = (data_70)*factor_70
print(data_70_corrected)
fits.writeto('ngc891_pacs_70_corrected.fits', data_70_corrected, header_70)


#Convolve:
convolved_image_70_to_350_new = convolve(data_70_corrected, gauss_70_to_350)
fits.writeto('70_convolved_to_350.fits', convolved_image_70_to_350_new, header_70)


#Regrid:
regridded_70 = mProject('ngc891_pacs_70_corrected.fits', '70_regridded.fits','target_70.txt')

regridded_70_to_350 = mProject('70_convolved_to_350.fits', '70_regridded_to_350.fits', 'target_350.txt')

data_regridded_70_to_350 = fits.getdata('70_regridded_to_350.fits')

header_regridded_70_to_350 = fits.getheader('70_regridded_to_350.fits')

data_regridded_70_to_350_scaled = data_regridded_70_to_350*pix_ratio_350_by_70

fits.writeto('70_regridded_to_350.fits', data_regridded_70_to_350_scaled, header_regridded_70_to_350, clobber=True)



#------------------------------------------------------------------------------
"""PACS 100 MICRON IMAGES - CONVOLVE/REGRID"""
#------------------------------------------------------------------------------


#Correct the data by multiplying with the conversion factor for the 100:
data_100_corrected = (data_100)*factor_100
print(data_100_corrected)
fits.writeto('ngc891_pacs_100_corrected.fits', data_100_corrected, header_100)


#Convolve:
convolved_image_100_to_350_new = convolve(data_100_corrected, gauss_100_to_350)
fits.writeto('100_convolved_to_350.fits', convolved_image_100_to_350_new, header_100)


#Regrid:
regridded_100 = mProject('ngc891_pacs_100_corrected.fits', '100_regridded.fits','target_100.txt')

regridded_100_to_350 = mProject('100_convolved_to_350.fits', '100_regridded_to_350.fits', 'target_350.txt')

data_regridded_100_to_350 = fits.getdata('100_regridded_to_350.fits')

header_regridded_100_to_350 = fits.getheader('100_regridded_to_350.fits')

data_regridded_100_to_350_scaled = data_regridded_100_to_350*pix_ratio_350_by_100

fits.writeto('100_regridded_to_350.fits', data_regridded_100_to_350_scaled, header_regridded_100_to_350, clobber=True)



#------------------------------------------------------------------------------
"""PACS 160 MICRON IMAGES - CONVOLVE/REGRID"""
#------------------------------------------------------------------------------


#Correct the data by multiplying with the conversion factor for the 70:
data_160_corrected = (data_160)*factor_160
print(data_160_corrected)
fits.writeto('ngc891_pacs_160_corrected.fits', data_160_corrected, header_160)


#Convolve:
convolved_image_160_to_350_new = convolve(data_160_corrected, gauss_160_to_350)
fits.writeto('160_convolved_to_350.fits', convolved_image_160_to_350_new, header_160)


#Regrid:
regridded_160 = mProject('ngc891_pacs_160_corrected.fits', '160_regridded.fits','target_160.txt')

regridded_160_to_350 = mProject('160_convolved_to_350.fits', '160_regridded_to_350.fits', 'target_350.txt')

data_regridded_160_to_350 = fits.getdata('160_regridded_to_350.fits')

header_regridded_160_to_350 = fits.getheader('160_regridded_to_350.fits')

data_regridded_160_to_350_scaled = data_regridded_160_to_350*pix_ratio_350_by_160

fits.writeto('160_regridded_to_350.fits', data_regridded_160_to_350_scaled, header_regridded_160_to_350, clobber=True)



#------------------------------------------------------------------------------
"""SPIRE 250 MICRON IMAGES - CONVOLVE/REGRID"""
#------------------------------------------------------------------------------


#Correct the data by multiplying with the conversion factor for the 70:
data_250_corrected = (data_250)*factor_250
print(data_250_corrected)
fits.writeto('ngc891_spire_250_corrected.fits', data_250_corrected, header_250)


#Convolve:
convolved_image_250_to_350_new = convolve(data_250_corrected, gauss_250_to_350)
fits.writeto('250_convolved_to_350.fits', convolved_image_250_to_350_new, header_250)


#Regrid:
regridded_250 = mProject('ngc891_spire_250_corrected.fits', '250_regridded.fits','target_250.txt')

regridded_250_to_350 = mProject('250_convolved_to_350.fits', '250_regridded_to_350.fits', 'target_350.txt')

data_regridded_250_to_350 = fits.getdata('250_regridded_to_350.fits')

header_regridded_250_to_350 = fits.getheader('250_regridded_to_350.fits')

data_regridded_250_to_350_scaled = data_regridded_250_to_350*pix_ratio_350_by_250

fits.writeto('250_regridded_to_350.fits', data_regridded_250_to_350_scaled, header_regridded_250_to_350, clobber=True)



#------------------------------------------------------------------------------
"""SPIRE 350 MICRON IMAGES - CONVOLVE/REGRID"""
#------------------------------------------------------------------------------


#Correct the data by multiplying with the conversion factor for the 350:
data_350_corrected = (data_350)*factor_350
print(data_350_corrected)
fits.writeto('ngc891_spire_350_corrected.fits', data_70_corrected, header_70)


#Regrid:
regridded_350 = mProject('ngc891_spire_350_corrected.fits', '350_regridded.fits', 'target350.txt.')



#------------------------------------------------------------------------------
"""REGRIDDED DATA"""
#------------------------------------------------------------------------------

data_1mm_regridded = fits.getdata('1mm_regridded.fits')
data_2mm_regridded = fits.getdata('2mm_regridded.fits')
data_70_regridded = fits.getdata('70_regridded.fits')
data_100_regridded = fits.getdata('100_regridded.fits')
data_160_regridded = fits.getdata('160_regridded.fits')
data_250_regridded = fits.getdata('250_regridded.fits')
data_350_regridded = fits.getdata('350_regridded.fits')


data_1mm_regridded_to_350 = fits.getdata('1mm_regridded_to_350.fits')
data_2mm_regridded_to_350 = fits.getdata('2mm_regridded_to_350.fits')
data_70_regridded_to_350 = fits.getdata('70_regridded_to_350.fits')
data_100_regridded_to_350 = fits.getdata('100_regridded_to_350.fits')
data_160_regridded_to_350 = fits.getdata('160_regridded_to_350.fits')
data_250_regridded_to_350 = fits.getdata('250_regridded_to_350.fits')



#------------------------------------------------------------------------------
"""GET THE RATIO IMAGES"""
#------------------------------------------------------------------------------


divided_data_1mm_over_350 = numpy.divide(data_1mm_regridded_to_350,data_350_regridded)
divided_data_2mm_over_350 = numpy.divide(data_2mm_regridded_to_350,data_350_regridded)
divided_data_70_over_350 = numpy.divide(data_70_regridded_to_350,data_350_regridded)
divided_data_100_over_350 = numpy.divide(data_100_regridded_to_350,data_350_regridded)
divided_data_160_over_350 = numpy.divide(data_160_regridded_to_350,data_350_regridded)
divided_data_250_over_350 = numpy.divide(data_250_regridded_to_350,data_350_regridded)


fits.writeto('1mm_by_350.fits', divided_data_1mm_over_350, header_350)
fits.writeto('2mm_by_350.fits', divided_data_2mm_over_350, header_350)
fits.writeto('70_by_350.fits', divided_data_70_over_350, header_350)
fits.writeto('100_by_350.fits', divided_data_100_over_350, header_350)
fits.writeto('160_by_350.fits', divided_data_160_over_350, header_350)
fits.writeto('250_by_350.fits', divided_data_250_over_350, header_350)


