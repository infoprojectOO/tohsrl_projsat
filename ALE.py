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
import methutil as methu
import sys
import graphics as gp
import numpy as np
import sgp4
import matplotlib.pyplot as mpl
import datetime
import methutil as methu
from matplotlib.patches import FancyArrowPatch
from itertools import product, combinations
import param  as cst
import subprocess
from satclass import *
from graphics import Painter
# del sys.modules['satclass']
# import satclass
# reload(satclass)

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

# mu_G = 398600.4418 * 10**9 # m³/s²
# R_G = 6378.1 * 10**3 # m 

# MM_air = 0.0289644 # kg/mol
# N_A = 6.023 * 10**23
# sigma_c = 1e-19 # m²


# Satellite orbital parameters
inputTLE = False
TLEfile = 'EGM96coefficients'
propagator = 'Classic'

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
    a_sat = R_G + 350 * 10**3 # semimajor axis - m
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

# # Projectile properties
# d_p = 0.01 #m
# rho_p = 8.96 * 1000 # kg/m³ copper

# Ejection parameters
eject_theta = 180. * m.pi/180. # (inverse) Pitching
eject_phi = 0. * m.pi/180. # (inverse) Yaw
eject_v_abs = 350 # m/s
eject_v_rel = 0.9

#Object creation

mechanics = Mechanics.initialise()
earth = Earth.create()
earth.setDate(date)
projectile = Projectile(cst.d_p,cst.Material.COPPER)
proj_spec = Projectile(cst.d_p,cst.Material.COPPER,Projectile.smooth)
proj_diff = Projectile(cst.d_p,cst.Material.COPPER,Projectile.coarse)
orbit = Orbit((a_sat,e_sat,i_sat,Om_sat,wp_sat))
schedule = Schedule()
timegen = TimeGenerators()
mechanics.set(earth,schedule)


if(inputTLE):
    nu0_sat = orbit.anomaly(M0_sat,'true')
ale_sat = Satellite(orbit,nu0_sat,R_G/10)

mechanics.add_animate(ale_sat)
mechanics.add_animate(earth)

schedule.plan('Orbiting phase', None, timelap, time = 0, duration = iter([timelap]))
schedule.plan('Projectile ejection', ale_sat.eject, projectile,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projectile,10**4))
# schedule.plan('Coarse Projectile ejection', ale_sat.eject, proj_diff,(eject_v_abs*0.95,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(proj_diff) )
# schedule.plan('Smooth Projectile ejection', ale_sat.eject, proj_spec,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(proj_spec) )

mechanics.start()

t = mechanics.timeline[1:]

r0_sat = orbit.getPos(ale_sat.nu_0)
r_sat = orbit.getPos(ale_sat.nu)
v_sat = orbit.getVel2(ale_sat.nu)

print('Satellite Starting position : ',r0_sat, )
print('Satellite Flight time : ',timelap*n/(2*m.pi))
print('Satellite Ending position : ',r_sat)
print('Satellite Ending speed : ',v_sat)


# print('Projectile ejection speed : ', np.linalg.norm(projectile.v_0))
# print('Projectile inbound speed : ', np.linalg.norm(projectile.v))
# print(projectile.v_0)

hypbox = mechanics.boxes[projectile]
ghostectile = hypbox.ghost
ghostbox = hypbox.ghostbox
dsmc = hypbox.dsmc
dsmc.abort_all()


# Plot orbit trajectory
ani = gp.plot(t,earth,ale_sat,[projectile,ghostectile],[mechanics.boxes[projectile],ghostbox],animation = True, globe = False, marks = (1,ghostbox.markindex))
#anighost = gp.plot(t,earth,ale_sat,[ghostectile],[mechanics.boxes[ghostectile]],animation = True, globe = False)
#boxes = [mechanics.boxes[proj_diff],mechanics.boxes[proj_spec]]
#ani = gp.plot(t,earth,ale_sat,[proj_diff,proj_spec],boxes,animation = False, globe = False)

def solve(projectile,atmosphere):
    global earth
    m = projectile.m
    S = projectile.S
    r0 = projectile.traj[0,:] # Starting position
    v0 = projectile.vel[0,:] # Starting velocity
    r_end = projectile.traj[-1,:] # Ending position

    lat, lon, name = Reference.get_refloc()
    transorbit = Orbit.retrieve(r0,v0)

    # --------------------- old not working junk !
    # x = v | z = r
    # x,z = symbols('x, z', real=True)
    # a,b = fit_atm(atmosphere)
    # eq = [ m*mu_G/(z*1000+R_G)**2  - 0.5*1*(a*exp(b*z))*S*x**2, 
    #        x - sqrt(mu_G/(R_G+z*1000)-29.56*10**6)]
    # sym = [x,z]
    # sol = nls(eq,sym)
    # ----------------------

    nu0, nu_end = transorbit.get_nu([np.linalg.norm(r0),50*1000+R_G])
    nu_scale = np.linspace(nu0,nu_end,1001)
    orbtraj = transorbit.getPos(nu_scale)
    orbvel = transorbit.getVel(nu_scale)
    h_scale = (np.linalg.norm(orbtraj,axis = 1)-R_G)*0.001
    rho_scale, h_scale = atmosphere.profile(h_scale)
    Fg = m*mu_G/(R_G+h_scale*1000)**2
    Fd = 0.5*1*S*rho_scale*np.linalg.norm(orbvel,axis=1)**2
    slope_scale = transorbit.get_slope(nu_scale)
    h_low = np.argmin(abs(h_scale-50))
    h_upp = np.argmin(abs(h_scale-120))
    resind = np.argmin(abs(Fd[h_upp:h_low]-(Fg[h_upp:h_low]*np.sin(-slope_scale[h_upp:h_low]))))+h_upp

    ht = h_scale[resind] # target altitude    


    intf = methu.integrate((Fd[:resind],(Fg*np.sin(-slope_scale))[:resind]),orbtraj[:resind])
    gain = mu_G/(R_G+h_scale[-1]) - mu_G/(R_G+ht)
    subtract = gain*(1-abs(intf[0])/abs(intf[1]))

    sol = (h_scale[resind], np.linalg.norm(orbvel[resind,:]), np.sqrt(np.linalg.norm(orbvel[0,:])**2+2*subtract))
    figforce = mpl.figure()
    lg,ld = mpl.plot(h_scale,Fg*np.sin(-slope_scale),h_scale,Fd)
    mpl.title('Drag force extrapolation vs. projectile weight projected on its trajectory')
    mpl.xlabel('Altitude (km)')
    mpl.ylabel('Force (N)')
    mpl.legend((lg,ld), ('Weight','Drag'), loc=1)

    figtraj = mpl.figure()
    axtraj = figtraj.add_subplot(111,projection='3d')
    renderer = Painter(figtraj,axtraj)
    #renderer.paint(transorbit)
    #renderer.paint(earth)
   
    axtraj.plot(orbtraj[:,0], orbtraj[:,1], orbtraj[:,2],'c-')
    axtraj.plot(projectile.traj[:,0],projectile.traj[:,1],projectile.traj[:,2],'r-')

    #-------------------------------- Debugging ---------------------------

    # Dh = h_scale[0]-h_scale[-1]
    # dr = (np.roll(orbtraj,1,axis=0)-orbtraj)[0:-1]
    # drnorm = np.linalg.norm(dr,axis=1)
    # Dh_approx = np.vdot(drnorm,np.sin(abs(slope_scale[:-1])))

    # nu_hist = transorbit.get_nu(np.linalg.norm(projectile.traj,axis=1))
    # r_hist = transorbit.getPos(nu_hist)
    # v_hist = transorbit.getVel(nu_hist)
    # h_hist = (np.linalg.norm(r_hist,axis=1)-R_G)*0.001
    # mpl.figure()
    # mpl.plot(h_hist,np.linalg.norm(projectile.vel,axis=1),'r-')
    # mpl.plot(h_hist,np.linalg.norm(v_hist,axis=1),'c-')

    #----------------------------------------------------------------------

    return sol, intf, subtract

#sol, intf, subtract = solve(projectile,earth.atm)
#---------------------------------------------------------
#debugger
print(dsmc.sentinel.is_alive())
sp_writer = dsmc.writer
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
