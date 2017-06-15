# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:44:14 2017

@author: ASiapan
"""
import math as m
from density import density
import scipy as sp
from importlib import reload
import sys

import numpy as np
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
from itertools import product, combinations
import param
import subprocess
import satclass
del sys.modules['satclass']
import satclass
reload(satclass)


global rel2abs, mu_G, R_G

def anomaly(theta, e, output = 'true'):
	res_anomaly = theta
	if(output=='true'):
		f = lambda x : x - e*m.sin(x) - theta
		ec_anomaly = sp.optimize.root(f,0,method='lm').x[0]
		res_anomaly = m.acos( (m.cos(ec_anomaly)-e) / (1 - e*m.cos(ec_anomaly)) )
	elif(output=='mean'):
		ec_anomaly = m.acos( (e+m.cos(theta)) / (1+e*m.cos(theta)) ) # eccentric anomaly
		res_anomaly = ec_anomaly - e*m.sin(ec_anomaly)
	return res_anomaly

def orb_pos(orb_param):
	(a,e,i,Om,w,theta) = orb_param
	x_plane = a*(m.cos(theta)-e)
	y_plane = a*(1-e**2)*m.sin(theta)
	r = a*(1-e**2)/(1 + e*m.cos(theta))
	return rel2abs.dot(np.array([x_plane,y_plane,0.]).transpose())

def orb_vel(orb_param):
	(a,e,i,Om,w,theta) = orb_param
	v_x = -m.sqrt(mu_G/(a*(1-e**2)))*m.sin(theta)
	v_y = m.sqrt(mu_G/(a*(1-e**2)))*(e+m.cos(theta))
	return rel2abs.dot(np.array([v_x,v_y,0.]).T)



def propel(r_ini, v_ini):
        r = r_ini
        v = v_ini
        traj = [r_ini]
        return traj,v


# -------------------------------------------------------------------------------
mpl.close("all")


mu_G = 398600 * 10**9 # m³/s²
R_G = 6378 * 10**3 # m 

# Satellite orbital parameters
a_sat = R_G + 400 * 10**3 # semimajor axis - m
e_sat = 0.4 # ellipticity
i_sat = 0. * m.pi/180 # inclination - rad
Om_sat = 90. * m.pi/180 # ascend node - rad
wp_sat = 90. * m.pi/180 # periapsis argument - rad
nu0_sat = 0. * m.pi/180 # starting anomaly  - rad

n = m.sqrt(mu_G/a_sat**3)

userinput = False

# Location and time coordinates
if(userinput):
	latlong = input('ALE Latitude , Longitude : ').split(',')
	latitude = float(latlong[0])
	longitude = float(latlong[1])
else:
	longitude = 0 * m.pi/180 
	latitude = 90 * m.pi/180

# orb_param = [a_sat,e_sat,i_sat,Om_sat,wp_sat,nu0_sat]

# Projectile properties
d_p = 0.002 #m
rho_p = 8.96 * 1000 # kg/m³ copper
C_D= 2.5
orbit = Orbit((a_sat,e_sat,i_sat,Om_sat,wp_sat))

ale_sat = Satellite(orbit,nu0_sat,R_G/10)

r0_sat = orbit.getPos(nu0_sat)

projectile = Projectile(d_p,rho_p)
# abs2rel = np.array([[m.cos(Om_sat)*m.cos(wp_sat)-m.sin(Om_sat)*m.cos(i_sat)*m.sin(wp_sat), m.sin(Om_sat)*m.cos(wp_sat)+m.cos(Om_sat)*m.cos(i_sat)*m.sin(wp_sat), m.sin(i_sat)*m.sin(wp_sat)],
# 						[-m.cos(Om_sat)*m.sin(wp_sat)-m.sin(Om_sat)*m.cos(i_sat)*m.cos(wp_sat), -m.sin(Om_sat)*m.sin(wp_sat)+m.cos(Om_sat)*m.cos(i_sat)*m.cos(wp_sat), m.sin(i_sat)*m.cos(wp_sat)],
# 						[m.sin(i_sat)*m.sin(Om_sat), -m.sin(i_sat)*m.cos(Om_sat), m.cos(i_sat)]])

# rel2abs = np.array([[m.cos(Om_sat)*m.cos(wp_sat)-m.sin(Om_sat)*m.cos(i_sat)*m.sin(wp_sat), -m.sin(wp_sat)*m.cos(Om_sat)-m.cos(wp_sat)*m.cos(i_sat)*m.sin(Om_sat), m.sin(i_sat)*m.sin(Om_sat)],
# 						[m.cos(wp_sat)*m.sin(Om_sat)+m.sin(wp_sat)*m.cos(i_sat)*m.cos(Om_sat), -m.sin(Om_sat)*m.sin(wp_sat)+m.cos(Om_sat)*m.cos(i_sat)*m.cos(wp_sat), -m.sin(i_sat)*m.cos(Om_sat)],
# 						[m.sin(i_sat)*m.sin(wp_sat), m.sin(i_sat)*m.cos(wp_sat), m.cos(i_sat)]])

# r0_sat = orb_pos(orb_param)

timelap = 0.5*2*m.pi / n # s


# Satellite equation of motion

ale_sat.clocktick(timelap)
r_sat = orbit.getPos(ale_sat.nu)
v_sat = orbit.getVel(ale_sat.nu)

# M_0 = anomaly(nu0_sat, e_sat,'mean')
# M = M_0 + n*(timelap)
# nu_sat = anomaly(M,e_sat,'true')
# orb_param[-1] = nu_sat

# r_sat = orb_pos(orb_param)
# v_sat = orb_vel(orb_param)

print('Satellite Starting position : ',r0_sat)
print('Satellite Flight time : ',timelap*n/(2*m.pi))
print('Satellite Ending position : ',r_sat)
ale_sat.eject(projectile,(200,0,m.pi))
# Plot Earth rim 

# phi = np.linspace(0, 2 * np.pi, 10)
# theta = np.linspace(0, np.pi, 10)
# xg = R_G * np.outer(np.cos(phi), np.sin(theta))
# yg = R_G * np.outer(np.sin(phi), np.sin(theta))
# zg = R_G * np.outer(np.ones(np.size(phi)), np.cos(theta))

# # Plot orbit trajectory

# dots = np.linspace(0,2*m.pi,100)
# x_ell = a_sat*(np.cos(dots)-e_sat)
# y_ell = a_sat*m.sqrt(1-e_sat**2)*np.sin(dots)
# z_ell = np.zeros(x_ell.size)
# orbit_traj = rel2abs.dot(np.array([x_ell,y_ell,z_ell]))

# # Plot satellite at final position
# sat_rim = R_G/10
# r = [-sat_rim*0.5, sat_rim*0.5]
# combi = np.array(list(product(r,r,r))).T
# combi1,combi2 = np.hsplit(combi,2)
# rows,cols = combi.shape
# combis = np.array(list(zip(combi1,combi2))).reshape(rows*2,int(cols/2))
# corners = np.ones((8,1)).dot(satpos)+np.array(list(product(r,r,r)))
# # x_sat = corners[:,0].reshape(8,1).dot(np.ones((1,8)))
# x_sat = np.ones(combis.shape)*r_sat[0]+combis
# y_sat = np.ones(combis.shape)*r_sat[1]+np.roll(combis,-2,axis=0)
# z_sat = np.ones(combis.shape)*r_sat[2]+np.roll(combis,-4,axis=0)
# # y_sat = corners[:,1].reshape(8,1).dot(np.ones((1,8)))
# # z_sat = corners[:,2].reshape(8,1).dot(np.ones((1,8)))

# # Plot satellite velocity at final position
# mag = 200
# xs,ys,zs = r_sat
# xv,yv,zv = mag*v_sat+r_sat
# varrow = Arrow3D([xs,xv],[ys,yv],[zs,zv], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")

renderer = Painter()

figearth = mpl.figure()
ax = figearth.add_subplot(111, projection='3d')
renderer.paint(orbit,ax)
renderer.paint(ale_sat,ax)
renderer.paint('earth',ax)
# ax.plot(orbit_traj[0,:],orbit_traj[1,:],orbit_traj[2,:],color = 'g')
# ax.plot_wireframe(xg, yg, zg)
# ax.plot_surface(x_sat,y_sat,z_sat,color = 'm')
# ax.add_artist(varrow)
# ax.set_xlabel('X : vernal point')
# ax.set_ylabel('Y : wake dir')
# ax.set_zlabel('Z : geo north dir')





