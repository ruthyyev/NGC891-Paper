from astropy.io import fits 
import numpy as np 
import os 
import matplotlib.pyplot as plt

#where to find the maps:
os.chdir('/Users/c1541417/Documents/891Paper/Maps/Convolved')


im_1mm = fits.getdata('1mm_convolved_to_350.fits')
im_2mm = fits.getdata('2mm_convolved_to_350.fits')
im_70 = fits.getdata('70_convolved_to_350.fits')
im_100 = fits.getdata('100_convolved_to_350.fits')
im_160 = fits.getdata('160_convolved_to_350.fits')
im_250 = fits.getdata('250_convolved_to_350.fits')
im_350 = fits.getdata('ngc891_nika2_350_corrected.fits')


point_1mm = im_1mm[176, 163]
point_2mm = im_2mm[176, 163]
point_70 = im_70[176, 163]
point_100 = im_100[176, 163]
point_160 = im_160[176, 163]
point_250 = im_250[176, 163]
point_350 = im_350[176, 163]

point_1mm_2 = im_1mm[177, 164]
point_2mm_2 = im_2mm[177, 164]
point_70_2 = im_70[177, 164]
point_100_2 = im_100[177, 164]
point_160_2 = im_160[177, 164]
point_250_2 = im_250[177, 164]
point_350_2 = im_350[177, 164]


"""
for i in range(163, 165):
    for j in range(175, 177):
        point_1mm = im_1mm[i, j]
        point_2mm = im_2mm[i, j]
"""
  
 
print im_1mm.sum() 


x1 = 1000
x2 = 2000
x70 = 70
x100 = 100
x160 = 160
x250 = 250
x350 = 350


plt.plot(x1, im_1mm.sum(), 'bo')
plt.plot(x2, im_2mm.sum(), 'bo')
#plt.plot(x70, point_70, 'bo')
#plt.plot(x100, point_100, 'bo')
#plt.plot(x160, point_160, 'bo')
##plt.plot(x250, point_250, 'bo')
#plt.plot(x350, point_350, 'bo')
plt.xscale('log')
plt.yscale('log')
plt.show()







