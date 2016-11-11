"""Description: Estimate the modified blackbody fit of a galaxy given certain 
parameters"""

import numpy as np
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy import optimize


#------------------------------------------------------------------------------

#Constants:

pi = const.pi 
c = const.c
h = const.h 
k = const.k 
dist = 9.6
mass = 1.6902675e38

#------------------------------------------------------------------------------

#Define the planck function:

def planck(v, temp):
    """
    Input: frequency values
    Returns: B_nu(v,temp), the Planck distribution function in SI units
    """
    numer = 2.*h*(v)**3/(c)**2
    denom = exp(h*c/(k*temp))-1.
    return numer/denom

#------------------------------------------------------------------------------

#Define the modified blackbody function: 

def mod(beta, dist, mass):
    """
    Input: beta, distance to galaxy, dust mass
    Returns: Flux density at given frequency
    """
    flux = (((v)**beta)/(dist)**2)*mass*planck
    return flux
    
#------------------------------------------------------------------------------
    