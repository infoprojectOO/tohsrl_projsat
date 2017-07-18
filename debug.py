# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:59:17 2017

@author: ASiapan
"""
import six
from density import density
import atmosphere.nrlmsise_00_header as atm
from atmosphere.nrlmsise_00 import gtd7
from six.moves import map, xrange, zip, reduce
import numpy as np
from matplotlib.axes import Axes, rcParams
import matplotlib.pyplot as mpl
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import colors as mcolors
from methutil import cube
from scipy.optimize import curve_fit
from satclass import *

#Input conversion
mpl.close("all")
#fig = mpl.figure()
#ax = fig.add_subplot(111,projection = '3d')
sat = ale_sat
proj = projectile
#board = ax
#cam = ale_sat.camera
#atmosphere = earth.atm
#atm = earth.atm
# orbit = Orbit((a_sat,e_sat,i_sat,Om_sat,wp_sat))

r_sat  = sat.traj[0,:]
r_sat = np.array([np.sqrt(2)/2*R_G,0,np.sqrt(2)/2*R_G])
r = np.linalg.norm(r_sat)
theta = m.acos(r_sat[2]/r)
phi = m.atan(r_sat[1]/r_sat[0])
as_legendre = sp.special.lpmv
" Find the gravitational potential at the desired point. "
grav_acc = 0.0 # Potential of the gravitational field at the stateVec location.
grav_pot = 0.0
for n in range(0, 10+1): # Go through all the desired orders and compute the geoid corrections to the sphere.
    term = 0. # Contribution to the potential from the current degree and all corresponding orders.
    for k in range(n+1): # Go through all the orders corresponding to the currently evaluated degree.
        norm = np.sqrt((2*n+1)*m.factorial(n-k)/m.factorial(n+k))
        term += norm * as_legendre(k,n,np.cos(theta)) * (earth.Cc[n][k]*np.cos(k*phi) + earth.Sc[n][k]*np.sin( k*phi ))
        print((n,k)," : ",norm *as_legendre(k,n,m.cos(theta)))
            
    grav_acc += -(n+1)*m.pow(earth.R/r, n+1)/r * term # Add the contribution from the current degree.
    grav_pot += m.pow(earth.R/r, n+1) * term
    
grav_acc *= earth.mu/earth.R # Final correction.
grav_pot *= earth.mu/earth.R

" Compute the acceleration due to the gravity potential at the given point. "
grav_accel = -(r_sat/r) * grav_pot/r
for deg in range(10):
    for order in range(deg+1):
        norm = np.sqrt((2*deg+1)/2*m.factorial(deg-order)/m.factorial(deg+order))
        print((deg,order),' : ', norm*sp.special.lpmv(order,deg,0))
        