# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:44:14 2017

@author: ASiapan
"""
import math as m
from density import density
import scipy as sp
import random as rand
from importlib import reload
import sys
import graphics as gp
import numpy as np
import sgp4
import matplotlib.pyplot as mpl
import datetime
import methutil as methu
from matplotlib.patches import FancyArrowPatch
from itertools import product, combinations
import param
import subprocess
from satclass import *
# del sys.modules['satclass']
# import satclass
# reload(satclass)

def projectileTimeGenerator(projectile, N_itmax = 10**4):
    dt = 10**(-3)*R_G/np.linalg.norm(projectile.v_0)
    Kn = 10
    n_it = 0
    alt = np.linalg.norm(projectile.r)-R_G

    while (Kn >= 1 and n_it <= N_itmax and alt>=0):

        alt = np.linalg.norm(projectile.r)-R_G
        n_volmean = density(alt*0.001)*N_A/(MM_air)
        lpm = 1/(m.sqrt(2)*n_volmean*sigma_c)
        Kn = lpm/projectile.d
        n_it += 1

        yield dt

def sphe2cart(r,lat,lon):
    x = r*np.cos(lat)*np.cos(lon)
    y = r*np.cos(lat)*np.sin(lon)
    z = r*np.sin(lat)
    return np.array([x,y,z])

def tle2params(filename):
    f = open(filename,'r')
    line1 = f.readline()
    line2 = f.readline()
    year, day, ballistic_star_coefficient = int(line1[18:20]), float(line1[20:32]), float(line1[33:43])
    inclination, right_ascension, eccentricity, argument_periapsis, mean_anomaly, mean_motion =\
        float(line2[8:16]), float(line2[17:25]), float("0." + line2[26:33]), float(line2[34:42]), float(line2[43:51]),\
        float(line2[52:63])
    semimajoraxis = m.pow(mu_G/mean_motion**2,1./3)
    data = {'a':semimajoraxis, 'e' : eccentricity, 'i' : inclination, 'Om' : right_ascension, 'wp' : argument_periapsis, 'M' : mean_anomaly, 
            'Y' : year, 'D' : day, 'B_' : ballistic_star_coefficient}
    return data

mu_G = 398600 * 10**9 # m³/s²
R_G = 6378 * 10**3 # m 

MM_air = 0.0289644 # kg/mol
N_A = 6.023 * 10**23
sigma_c = 1e-19 # m²


# Satellite orbital parameters
inputTLE = False
TLEfile = 'EGM96coefficients'

if(inputTLE):
    params = tle2params(TLEfile)
    a_sat = params['a']
    e_sat = params['e']
    i_sat = params['i']
    Om_sat = params['Om']
    wp_sat = params['wp']
    M0_sat = params['M']
    year_sat = params['Y']
    day_sat = params['D']
    elapsed = ((year_sat-2000)*365.25 + day_sat)*86400
else:
    a_sat = R_G + 400 * 10**3 # semimajor axis - m
    e_sat = 0.0 # ellipticity
    i_sat = 0. * m.pi/180 # inclination - rad
    Om_sat = 0. * m.pi/180 # ascend node - rad
    wp_sat = 0. * m.pi/180 # periapsis argument - rad
    nu0_sat = 0. * m.pi/180 # starting anomaly  - rad
    date = datetime.datetime(2000, 1 , 1, 12, 0, 0, 0)
    dateref = datetime.datetime(2000,1,1,12)
    year_sat = date.year # reference year - years
    day_sat = date.day # day of the year - days
    elapsed = (date - dateref).total_seconds()

n = m.sqrt(mu_G/a_sat**3)

timelap = 0.25*2*m.pi / n # s

# Projectile properties
d_p = 0.002 #m
rho_p = 8.96 * 1000 # kg/m³ copper
C_D= 2.5

# Ejection parameters
eject_theta = 180. * m.pi/180. # (inverse) Pitching
eject_phi = 0. * m.pi/180. # (inverse) Yaw
eject_v_abs = 200 # m/s
eject_v_rel = 0.9

#Object creation

mechanics = Mechanics.initialise()
earth = Earth()
earth.setDate(date)
projectile = Projectile(d_p,rho_p)
orbit = Orbit((a_sat,e_sat,i_sat,Om_sat,wp_sat))
schedule = Schedule()
mechanics.set(earth,schedule)


if(inputTLE):
    nu0_sat = orbit.anomaly(M0_sat,'true')
ale_sat = Satellite(orbit,nu0_sat,R_G/10)

mechanics.add_animate(ale_sat)
mechanics.add_animate(earth)

schedule.plan('Orbiting phase', None, timelap, time = 0, duration = iter([timelap]))
schedule.plan('Projectile ejection', ale_sat.eject, projectile,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = projectileTimeGenerator(projectile) )

mechanics.start()

t = mechanics.timeline[1:]

r0_sat = orbit.getPos(ale_sat.nu_0)
# ale_sat.clocktick(timelap)
r_sat = orbit.getPos(ale_sat.nu)
v_sat = orbit.getVel2(ale_sat.nu)

print('Satellite Starting position : ',r0_sat, )
print('Satellite Flight time : ',timelap*n/(2*m.pi))
print('Satellite Ending position : ',r_sat)
print('Satellite Ending speed : ',v_sat)


print('Projectile ejection speed : ', np.linalg.norm(projectile.v_0))
print('Projectile inbound speed : ', np.linalg.norm(projectile.v))
print(projectile.v_0)


# Plot orbit trajectory
ani = gp.plot(t,earth,ale_sat,projectile,animation = True, globe = False)



#---------------------------------------------------------
#debugger

v_0 = projectile.v_0
nu = ale_sat.nu
M = ale_sat.M
vel2rel = orbit.getVel2Rel(nu)
dvx = eject_v_abs*(m.cos(eject_theta))
dvy = eject_v_abs*(m.sin(eject_theta)*m.cos(eject_phi))
dvz = eject_v_abs*(m.sin(eject_theta)*m.sin(eject_phi))
dv = np.array([dvx,dvy,dvz])
rel2abs = orbit.rel2abs
vel2rel = orbit.getVel2Rel(nu)

r = a_sat*(1-e_sat**2)/(1+e_sat*m.cos(nu))
v = m.sqrt(mu_G*(2/r-1/a_sat))
rnorm = np.linalg.norm(r_sat)
vnorm = np.linalg.norm(v_sat)
v_radmeth = orbit.getVel(nu)
v_tanmeth = orbit.getVel2(nu)

v_x = -m.sqrt(mu_G/(a_sat*(1-e_sat**2)))*m.sin(nu)
v_y = m.sqrt(mu_G/(a_sat*(1-e_sat**2)))*(e_sat+m.cos(nu))
v_r = m.sqrt(mu_G/(a_sat*(1-e_sat**2)))*e_sat*m.sin(nu)
v_t = m.sqrt(mu_G/(a_sat*(1-e_sat**2)))*(1+e_sat*m.cos(nu))

ref = np.array([0.,0.,1.])
pointing = r_sat/np.linalg.norm(r_sat)
rot_vec = np.cross(ref,pointing)
attitude = Quaternion(rot_vec/np.linalg.norm(rot_vec),m.asin(np.linalg.norm(rot_vec)))

xd, yd, zd = np.array(ale_sat.attitude.to_matrix()).dot(np.array([0,0,ale_sat.width]))+r_sat