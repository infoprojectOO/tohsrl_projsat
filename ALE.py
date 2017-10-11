# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:44:14 2017

@author: ASiapan

Main file for simulating projectile atmospheric reentry along with satellite orbiting and tracking.

"""
import math as m
import scipy as sp
import random as rand
from importlib import reload
import methutil as methu
import sys
import graphics as gp
import numpy as np
#import sgp4
import traceback
import matplotlib.pyplot as mpl
import datetime
import methutil as methu
import param  as cst
from satclass import *
from physenv import *
from ParametersEditor import ParamInput
from ds2v2py import manSimIter
# del sys.modules['satclass']
# import satclass
# reload(satclass)

def tle2params(filename):
    """ Function reads data from a file and stores it in a data structure object
    """
    f = open(filename,'r')
    line1 = f.readline()
    line2 = f.readline()
    year, day, ballistic_star_coefficient = int(line1[18:20]), float(line1[20:32]), float(line1[33:43])
    inclination, right_ascension, eccentricity, argument_periapsis, mean_anomaly, mean_motion =\
        float(line2[8:16]), float(line2[17:25]), float("0." + line2[26:33]), float(line2[34:42]), float(line2[43:51]),\
        float(line2[52:63])
    semimajoraxis = m.pow(cst.mu_G/mean_motion**2,1./3)
    data = {'a':semimajoraxis, 'e' : eccentricity, 'i' : inclination, 'Om' : right_ascension, 'wp' : argument_periapsis, 'M' : mean_anomaly, 
            'Y' : year, 'D' : day, 'B_' : ballistic_star_coefficient}
    return data

# cst.mu_G = 398600.4418 * 10**9 # m³/s²
# R_G = 6378.1 * 10**3 # m 

# MM_air = 0.0289644 # kg/mol
# N_A = 6.023 * 10**23
# sigma_c = 1e-19 # m²


# Satellite orbital parameters
inputTLE = False
TLEfile = 'EGM96coefficients'
propagator = 'Classic' # choose between classical keplerian orbits and gravity + drag corrected propagator for the satellite (to be implemented yet ! --- very straightforward though )

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
    a_sat = cst.R_G + 400 * 10**3 # semimajor axis - m
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

n = m.sqrt(cst.mu_G/a_sat**3)

timelap = 0.25*2*m.pi / n # s

# Ejection parameters
eject_theta = 180. * m.pi/180. # (inverse) Pitching
eject_phi = 0. * m.pi/180. # (inverse) Yaw
eject_v_abs = 350 # m/s
eject_v_rel = 0.9 # relative velocity fraction from satellite (not really used)

#Object creation

    #Physical environment
earth = Earth.create()
earth.setDate(date)
orbit = Orbit((a_sat,e_sat,i_sat,Om_sat,wp_sat))
if(inputTLE):
    nu0_sat = orbit.anomaly(M0_sat,'true')
    # Satellite
ale_sat = Satellite(orbit,nu0_sat,cst.R_G/10)
    # Time-management
mechanics = Mechanics.initialise()
schedule = Schedule()
timegen = TimeGenerators()
mechanics.set(earth,schedule)
mechanics.add_animate(ale_sat)
mechanics.add_animate(earth)

    # Projectile(s)
projectile = Projectile(cst.d_p,cst.Material.COPPER, ablation = 0); #projectile.add_spec(DynamicBox.OVERRIDE_DRAG, 2)
# proj_spec = Projectile(cst.d_p,cst.Material.COPPER,Projectile.smooth)
# proj_diff = Projectile(cst.d_p,cst.Material.COPPER,Projectile.coarse)
#proj0 = Projectile(cst.d_p,cst.Material.COPPER); proj0.add_spec(DynamicBox.OVERRIDE_AERODOM, AeroDom.Free_Molecular_Flow)#proj0.add_spec(DynamicBox.OVERRIDE_DRAG, 2)
proj1 = Projectile(cst.d_p,cst.Material.COPPER); proj1.add_spec(DynamicBox.OVERRIDE_DRAG, 1)
proj2 = Projectile(cst.d_p,cst.Material.COPPER); proj2.add_spec(DynamicBox.OVERRIDE_DRAG, 2)
# projectiles = [projectile,proj1,proj2]
projCu = Projectile(cst.d_p,cst.Material.COPPER,finition = 0.8, ablation = 0.01, id = 'Cu')
projZn = Projectile(cst.d_p,cst.Material.ZINC,finition = 0.8, ablation = 0.01, id = 'Zn')
projFe = Projectile(cst.d_p,cst.Material.IRON,finition = 0.8, ablation = 0.01, id = 'Fe')
projTi = Projectile(cst.d_p,cst.Material.TITANIUM,finition = 0.8, ablation = 0.01, id = 'Ti')
# projs = [projCu,projZn,projFe,projTi]
projperm = Projectile(cst.d_p,cst.Material.COPPER,finition = 0.8, ablation = 0.01)
projmelt = Projectile(cst.d_p,cst.Material.COPPER,finition = 0.8, ablation = 0)

projs = [projFe] # List of all projectiles used (launched)

schedule.plan('Orbiting phase', None, timelap, time = 0, duration = iter([timelap]))
# projectile = Projectile(cst.d_p,cst.Material.COPPER,finition = Projectile.coarse,ablation=0)
# projectile.add_spec(DynamicBox.OVERRIDE_DRAG, 2)
# schedule.plan('Projectile ejection Cu', mechanics.eject, ale_sat, projCu,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projCu) )
# schedule.plan('Projectile ejection Zn', mechanics.eject, ale_sat, projZn,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projZn) )
# schedule.plan('Projectile ejection Fe', mechanics.eject, ale_sat, projFe,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projFe) )
# schedule.plan('Projectile ejection Ti', mechanics.eject, ale_sat, projTi,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projTi) )
schedule.plan('Projectile ejection perm', mechanics.eject,ale_sat, projperm,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projperm) )
schedule.plan('Projectile ejection melt', mechanics.eject,ale_sat, projmelt,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projmelt), id = 'melt' )

# schedule.plan('Projectile ejection', mechanics.eject, ale_sat, projectile,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(projectile))
#schedule.plan('Projectile ejection 0', ale_sat.eject, proj0,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(proj0) )
# schedule.plan('Projectile ejection 2', ale_sat.eject, proj2,(eject_v_abs,eject_theta,eject_phi), time = timelap, duration = timegen.projectileTimeGenerator(proj2) )

try : # start the universe clock
    mechanics.start()
except:
    traceback.print_exc()
finally: # Retrieve all data computed so far and abort all running threads
    t = mechanics.timeline[1:]
    for p in mechanics.boxes:
        if p.__class__.__name__ == 'Projectile':
            hypbox = mechanics.boxes[p]
            ghostectile = hypbox.ghost
            ghostbox = hypbox.ghostbox
            dsmc = hypbox.dsmc
            db = dsmc.database
            dsmc.abort_all()
            dsmc.save_data()

r0_sat = orbit.getPos(ale_sat.nu_0)
r_sat = orbit.getPos(ale_sat.nu)
v_sat = orbit.getVel2(ale_sat.nu)

print('Satellite Starting position : ',r0_sat)
print('Satellite Ending position : ',r_sat)
print('Satellite Ending speed : ',v_sat)


# Plot orbit trajectory

# boxes = [mechanics.boxes[p] for p in projs]+[ghostbox]
# palette = gp.palette(0,len(projs),natural = True)
# ani, axdic = gp.plot(t,earth,ale_sat,projs+[ghostectile],boxes,animation = False, save = False, globe = False, dsmcdata = dsmc.database.values(),  marks = (1,ghostbox.markindex))
#ani, axdic = gp.plot(t,earth,ale_sat,[projectile],[mechanics.boxes[projectile]],animation = True, save = False, globe = False)
#ani, axdic = gp.plot(t,earth,ale_sat,[projectile,ghostectile],[mechanics.boxes[projectile],ghostbox],animation = True, save = False, globe = False, marks = (1,ghostbox.markindex), dsmcdata = db.values())
ani, axdic = gp.plot(t,earth,ale_sat,[projperm,projmelt],[mechanics.boxes[projperm],mechanics.boxes[projmelt]],animation = False, save = False, globe = False)
#ani = gp.plot(t,earth,ale_sat,[proj0],[mechanics.boxes[proj0]],animation = True)
# ani = gp.plot(t,earth,ale_sat,[proj1,proj2],[mechanics.boxes[proj1],mechanics.boxes[proj2]],animation = True, globe = False)
#anighost = gp.plot(t,earth,ale_sat,[ghostectile],[mechanics.boxes[ghostectile]],animation = True, globe = False)
#boxes = [mechanics.boxes[proj_diff],mechanics.boxes[proj_spec]]
#ani = gp.plot(t,earth,ale_sat,[proj_diff,proj_spec],boxes,animation = False, globe = False)
# gp.compare(ale_sat,[proj1,proj2],[mechanics.boxes[proj1],mechanics.boxes[proj2]])


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
v = m.sqrt(cst.mu_G*(2/r-1/a_sat))
rnorm = np.linalg.norm(r_sat)
vnorm = np.linalg.norm(v_sat)
v_radmeth = orbit.getVel(nu)
v_tanmeth = orbit.getVel2(nu)

v_x = -m.sqrt(cst.mu_G/(a_sat*(1-e_sat**2)))*m.sin(nu)
v_y = m.sqrt(cst.mu_G/(a_sat*(1-e_sat**2)))*(e_sat+m.cos(nu))
v_r = m.sqrt(cst.mu_G/(a_sat*(1-e_sat**2)))*e_sat*m.sin(nu)
v_t = m.sqrt(cst.mu_G/(a_sat*(1-e_sat**2)))*(1+e_sat*m.cos(nu))

ref = np.array([0.,0.,1.])
pointing = r_sat/np.linalg.norm(r_sat)
rot_vec = np.cross(ref,pointing)
attitude = Quaternion(rot_vec/np.linalg.norm(rot_vec),m.asin(np.linalg.norm(rot_vec)))

xd, yd, zd = np.array(ale_sat.attitude.to_matrix()).dot(np.array([0,0,ale_sat.width]))+r_sat


# path = r'D:\Program\DSMC\DS2V\DS2VD.dat'
# lat, lon = earth.pos2coord(projectile.r)
# air = earth.atm.at(0,np.linalg.norm(projectile.r)-R_G,lat, lon)
# print(air.p)

# print('Starting contdition:\nDensity = ' + str(air.nv) + ',\t Temp = ' + str(air.T) + ',\t Vx = ' + str(np.linalg.norm(projectile.v)) + '\nat alt = ' + str((np.linalg.norm(projectile.r)-R_G)/1000))
# ParamInput(air.nv,air.T,air.T,np.linalg.norm(projectile.v),path)

# print('\nNow please open DS2V for simulation')
# check = ''
# cont = 'y'
# i_it = 0
# while (check != 'r'):
#     check = input('Ready? (r): ')
# while (cont == 'y'):
#     Cd = input('Please retrieve Drag Coef from DS2V\n Cd = ')
#     rho, nv, T, v_pnorm, alt, projectile = manSimIter(projectile,earth,float(Cd))
#     ParamInput(nv,T,T,v_pnorm,path)
#     print('Density = ' + str(nv) + ',\t Temp = ' + str(T) + ',\t Vx = ' + str(v_pnorm) + '\nat alt = ' + str(alt/1000))
#     print('Surface:' + str(projectile.S))
#     cont = input('Continue? (y/n): ')

