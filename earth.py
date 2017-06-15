# -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:06:42 2017

@author: ASiapan
"""
import PIL
import matplotlib.pyplot as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Spatial mechanics parameters

mu_G = 398600
R_G = 6378
r_p = 100+R_G
r_a = 400+R_G
a_p = 250+R_G

# Projectile properties
d_p = 0.002 # m
S_p = m.pi*d_p**2
C_D= 2.5
v_p = m.sqrt(2*mu_G*((1/r_p)-1/(2*a_p)))*1000 # m/s

print(v_p)

# Air properties
alt_ref = 71 #km
dst_ref = 6.4*10**(-5) #kg/m^3
T_ref = 214.65 #K
c_p = 1000 #J/kg K
sigma = 5.67*10**(-8) #W/m²K^4
k_a = 0.05 #W/mK
M = 0.0289644 #kg/mol
R = 8.3144 #J/Kmol
k_B = 1.38*10**(-23) # J/K
R_M = R/M
g0 = 9.81 #m/s^2
gamma = 1.4
print(dst_ref*m.e**(-(100-alt_ref)*1000*g0*M/(T_ref*R)))

T_a = T_ref + v_p**2/(2*c_p)
T_0 = 8000 #K

c_a = m.sqrt(gamma*R_M*223.15)
rho = density(100.0)
print("density at 100 km : %f" % rho)

v_therm = m.sqrt(gamma*R_M * T_ref)


F_drag = 0.5*C_D*rho*v_p**2 #N
P_rad = F_drag*v_p #W

t_c = (c_p*rho)/(3*T_0**3*sigma)

k_a/(c_p*rho*0.1**2)
sigma*(8000-T_ref)**3

print("Average power radiated at peak emission : %f" % P_rad)

# -------------------------------------------------------------------------------

# load bluemarble with PIL
bm = PIL.Image.open('cylmapearth.jpg')
# it's big, so I'll rescale it, convert to array, and divide by 256 to get RGB values that matplotlib accept 
bm = np.array(bm)/256.

# coordinates of the image - don't know if this is entirely accurate, but probably close
lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 

# repeat code from one of the examples linked to in the question, except for specifying facecolors:
fig = mpl.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.outer(np.cos(lons), np.cos(lats)).T
y = np.outer(np.sin(lons), np.cos(lats)).T
z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors = bm)

mpl.show()
