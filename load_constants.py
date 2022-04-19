#*****************************************************************************#
#                   Constants for Sim                                         #
#*****************************************************************************#

#---
#  Physical Constants and Earth Shaping Constants
#---
import numpy as np
mu_E =3.986004418e5  #earth's gravitational parameters 
r_E =6.3781e3 #earth's radius, km
J2 = 0
J3 = 0 
J4 = 0
J5 = 0
# spherical harmonics 

#---
# Simulation options
#---
delta_t = 0.01 #step size
save_data = 10.0 # when to save to  a data frame
