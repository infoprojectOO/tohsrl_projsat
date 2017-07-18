# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:20:58 2017

@author: ASiapan
"""

# Class containing moving objects and their properties

import math as m
import numpy as np
import scipy as sp
import csv
import atmosphere.nrlmsise_00_header as atm
from atmosphere.nrlmsise_00 import gtd7
import datetime as dt
from density import density
from mathutils import Quaternion, Vector
import methutil as methu
global mu_G, R_G
mu_G = 398600.4418 * 10**9 # m³/s²
R_G = 6378.1 * 10**3 # m 

class Movody:
    '''Generic class for moving bodies = (movodies) in space'''

    def __init__(self, mass, surface, position0, velocity0):
        self.m = mass
        self.S = surface
        self.r = position0.copy()
        self.v = velocity0.copy()
        self.traj = np.array([position0]) 
        self.vel = np.array([velocity0])

    def clocktick(self,dv,dt):
        self.r += self.v*dt + dv*dt*0.5 #Euler forward integration
        self.v += dv
        self.traj = np.vstack((self.traj,self.r))
        self.vel = np.vstack((self.vel,self.v))

class Satellite(Movody):
    ''' Class Satellite captures the properties and abilities of the satellite
    It moves along an orbit, ejects and tracks a projectile and possesses a camera'''
    
    def __init__(self, orbit, nu_ini, width):
        super().__init__(5,width,orbit.getPos(nu_ini),orbit.getVel(nu_ini))
        self.orbit = orbit
        self.tracking = False # no projectile to track yet
        self.nu_0 = self.nu = nu_ini # real anomaly
        self.M_0 = self.M = orbit.anomaly(nu_ini,'mean') # mean anomaly
        self.traj = np.array([orbit.getPos(nu_ini)]) 
        self.speed = np.array([orbit.getVel(nu_ini)])
        self.width = self.d = width
        self.camera = Camera(width*0.8)
        self.att0 = Quaternion([0.,1.,0.],m.pi*0.5)
        self.setattitude()
        self.attitudes = [self.attitude.copy()]
        

    def clocktick(self, dt):
        self.M += self.orbit.n*dt
        self.M = self.M%(2*m.pi)
        self.nu = self.orbit.anomaly(self.M,'true')
        self.setattitude()
        if(self.tracking):
            self.track()


    def eject(self, projectile, dv):
        (norm,theta,phi) = dv
        dvx = norm*(m.cos(theta)) # x points towards the satellite's direction
        dvy = norm*(m.sin(theta)*m.cos(phi)) # y points inwards in the orbital plane
        dvz = norm*(m.sin(theta)*m.sin(phi))  # z completes the triorthon
        r = self.orbit.getPos(self.nu)
        v = self.orbit.getVel(self.nu)+self.orbit.rel2earth.dot(self.orbit.getVel2Rel(self.nu).dot(np.array([dvx,dvy,dvz]).T)).T
        projectile.propel(r,v)
        self.setTarget(projectile)
        Mechanics.getMech().add_animate(projectile)
        return projectile

    def setTarget(self, projectile):
        self.obj = projectile
        self.tracking = True
        r = self.orbit.getPos(self.nu)
        self.linesight = projectile.r-r
        if(np.linalg.norm(self.linesight)!=0.):
            self.linesight = self.linesight/np.linalg.norm(self.linesight)
        self.lostraj = np.array([self.linesight.copy()])
        self.traj = np.array([r.copy()])
        #self.setattitude(-self.linesight)

    def track(self):
        r = self.orbit.getPos(self.nu)
        v = self.orbit.getVel(self.nu)
        self.linesight = self.obj.r-r
        self.linesight = self.linesight/np.linalg.norm(self.linesight)
        self.lostraj = np.vstack((self.lostraj,self.linesight))
        self.traj = np.vstack((self.traj,r))
        self.speed = np.vstack((self.speed,v))
        self.attitudes = self.attitudes+[self.attitude]
        #self.setattitude(-self.linesight)


    def setattitude(self):
        #ref = Vector([1.,0.,0.])
        self.attitude = self.att0.copy()
        Om_rot = Quaternion([0.,0.,1.],self.orbit.Om)
        i_rot = Quaternion([1.,0.,0.],self.orbit.i)
        i_rot.rotate(Om_rot)
        wnu_rot = Quaternion([0.,0.,1.],self.orbit.w+self.nu)
        wnu_rot.rotate(i_rot)        
        self.attitude.rotate(wnu_rot)
        # pointing = np.array(orbit.getPos(self.nu)/np.linalg.norm(orbit.getPos(self.nu)))
        # rot_vec = np.cross(ref,pointing)
        # rot_norm = np.linalg.norm(rot_vec)
        # self.attitude = Quaternion(rot_vec/rot_norm,m.asin(rot_norm))

class Camera:

    def __init__(self,width):
        self.width = width

class Projectile(Movody):
    ''' Passive object representing a solid sphere and recording its kinetic state trajectory '''

    def __init__(self, diameter, density):
        self.d = diameter
        self.S = 0.25*m.pi*diameter**2
        self.m = m.pi/6.*diameter**3*density

    def propel(self, r_ini, v_ini):
        self.r_0 = r_ini
        self.r = r_ini.copy()
        self.v_0 = v_ini
        self.v = v_ini.copy()
        self.traj = np.array([r_ini])
        self.vel = np.array([v_ini])

    def clocktick(self,dt,dv):
        self.r += self.v*dt + dv*dt*0.5
        self.v += dv
        self.traj = np.vstack((self.traj,self.r))
        self.vel = np.vstack((self.vel,self.v))