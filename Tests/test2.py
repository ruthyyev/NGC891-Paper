from astropy.io import fits 
import os

os.chdir('/Users/c1541417/Documents/891Paper/Maps')

data, header = fits.getdata('ngc891_spire_350_1.fits', header = True)

header['XTENSION'] = 'IMAGE'  
header['BITPIX'] = -64                                                                                                                                 
header['PCOUNT'] = 0                         
header['GCOUNT'] = 1                                      
header['LONGSTRN'] = 'OGIP 1.0'                                          
header['META_0']  = 327                                         
header['META_1']  = 351                                          
header['CRPIX1']  = 163.0                                                               
header['CRPIX2'] = 179.0                                                             
header['CRVAL1']  =    35.64035037478851  
header['CRVAL2']  =    42.35599876715473  
header['CDELT1']  =   -0.002222222222222        
header['CDELT2']  =    0.002222222222222       
header['CTYPE1']  = 'RA---TAN' 
header['CTYPE2']  = 'DEC--TAN' 
header['EQUINOX'] = 2000.0                       
header['CROTA2'] = 0.0                      

fits.writeto('ngc891_spire_350_new2.fits', data, header, clobber = True)
