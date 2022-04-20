#*****************************************************************************#
#                   Constants for Sim                                         #
#*****************************************************************************#

#---
#  Physical Constants and Earth Shaping Constants
#---
import numpy as np
mu_E =3.986004418e5  #earth's gravitational parameters 
r_E =6.3781e3 #earth's radius, km
J2 = 1.0826269e-3 # From Vallado probably outdated
J3 = -2.5323e-6
J4 = -1.6204e-6
J5 = 5.40676e-7
# spherical harmonics 

#---
# Simulation options
#---
delta_t = 0.01 #step size
save_data = 10.0 # when to save to  a data frame
