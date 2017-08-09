# -*- coding: utf-8 -*-
"""
Created on Thu May 25 11:42:35 2017

@author: ASiapan
"""

import math as m
import param as cst
import numpy as np
import scipy as sp
import csv
import time
import os, pathlib, shlex
import subprocess
import atmosphere.nrlmsise_00_header as atm
from atmosphere.nrlmsise_00 import gtd7
import datetime as dt
from density import density
from mathutils import Quaternion, Vector
import methutil as methu
from threading import Thread
import threading
import dsmc

# global cst.mu_G, R_G
# cst.mu_G = 398600.4418 * 10**9 # m³/s²
# R_G = 6378.1 * 10**3 # m 

# k_B = 1.38064852*10**(-23) #J/K
# MM_air = 0.0289644 # kg/mol
# N_A = 6.023 * 10**23
# sigma_c = 1e-19 # m²7
# sigma_r = 5.670373 * 10**(-8) # W/m²K^4
# k_c = 401 # W/(m K)

class Singleton(type):
    ''' Implementation of the Singleton design pattern on Python'''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        # print('called by :',cls)
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        # print(Singleton._instances)
        return cls._instances[cls]

class Earth(metaclass=Singleton):
    ''' Earth class encompasses the natural elements proper to earth :
    * Time frame(s) : Sidereal and Solar
    * Geometry,
    * Gravity field,
    * Atmosphere
    '''
    R = cst.R_G
    mu = cst.mu_G

    def __init__(self):
        self.sidereal = 86164.1 # sidereal day duration
        self.solarday = 86400 # solar day duration
        self.w_rot = 2*m.pi/self.sidereal
        self.tilt = 23.44*m.pi/180
        self.tilt_rot = Quaternion([1,0,0],-self.tilt)
        self.axis = self.tilt_rot*Vector((0,0,1))
        self.rotang = 0.
        self.clock = 0.
        self.sunclock = 0.
        self.phi = [self.rotang]
        self.Cc, self.Sc = self.unpackgravitymodel()
        self.MAX_DEG = 10

    @classmethod
    def create(cls):
        ''' Method for replacing the old mechanics object by a new one, must by imperatively called first if program in a script'''
        cls._instances.pop(cls,'empty')
        return Earth()

    def setDate(self,date):
        self.clock = (date - Reference.get_refdate()).total_seconds()%self.sidereal
        self.sunclock = (date - Reference.get_refdate()).total_seconds()%self.solarday
        self.rotang = (self.clock*self.w_rot)%2*m.pi
        self.phi = [self.rotang]
        self.atm = Atmosphere(date.year, date.day, model = 'advanced')

    def clocktick(self,dt):
        self.clock += dt
        self.sunclock += dt
        self.sunclock = self.sunclock%self.solarday
        self.rotang = ((self.clock%self.sidereal)*self.w_rot)%2*m.pi
        self.phi.append(self.rotang)

    def pos2coord(self,r):
        x_rel, y_rel, z_rel = self.abs2rel(r)
        #rnorm = np.linalg.norm(r)

        lon = m.atan2(y_rel,x_rel)
        lat = m.atan2(z_rel,m.sqrt(y_rel**2+x_rel**2))

        return lat, lon

    def abs2rel(self,r):
        rel_rot = -self.tilt_rot
        rel_rot.rotate(Quaternion(self.axis,-self.rotang))
        g_rot = np.array(rel_rot.to_matrix())
        return g_rot.dot(r)

    def unpackgravitymodel(self):
        """ Read the EGM96 gravitational field model coefficients from EGM96coefficients
        file and parse them to be used with computeGravitationalPotential functions.
        
        Returns
        ----------
        2-tuple of the C and S coefficients of EGM96 model. They are stored in dictionaries
            of list. The keys are degrees of the potential expansion and the values()
            of the list entries are the coefficients corresponding to the orders for
            the expansion to a given degree.
        
        Reference
        ----------
        EGM96 coefficients have been downloaded from:
            ftp://cddis.gsfc.nasa.gov/pub/egm96/general_info/readme.egm96
        """
        " Read the coefficients. "
        Ccoeffs = {0:[1],1:[0,0]}; Scoeffs ={0:[0],1:[0,0]};
        degrees = [];
        with open("EGM96coefficients", "r") as egm96file:
            reader = csv.reader(egm96file, delimiter=" ")
            olddeg = 1
            for row in reader:
                newdeg = int(row[1])
                if(newdeg>olddeg):
                    degrees.append( int(row[1]) ) # convert strings into numbers
                    Ccoeffs[newdeg] = [float(row[3])]
                    Scoeffs[newdeg] = [float(row[4])]
                else:
                    Ccoeffs[newdeg].append(float(row[3]))
                    Scoeffs[newdeg].append(float(row[4]))
                olddeg = newdeg

        return Ccoeffs, Scoeffs

    def gravity(self,r_sat):
        """ Calculate the acceleration due to gravtiy acting on the satellite at
        a given state (3 positions and 3 velocities). Ignore satellite's mass,
        i.e. use a restricted two-body problem.
        Arguments
        ----------
        numpy.ndarray of shape (1,6) with three Cartesian positions and three
            velocities in an inertial reference frame in metres and metres per
                second, respectively.
        epoch - float corresponding to the epoch at which the rate of change is
            to be computed.
        Returns
        ----------
        numpy.ndarray of shape (1,3) with three Cartesian components of the
            acceleration in m/s2 given in an inertial reference frame.
        """
        " Compute geocentric latitude and longitude. "
        r_rel = self.abs2rel(r_sat)
        r = np.linalg.norm(r_rel)
        theta = m.acos(r_rel[2]/r)
        phi = m.atan(r_rel[1]/r_rel[0])
        as_legendre = sp.special.lpmv
        " Find the gravitational potential at the desired point. "
        grav_acc = 0.0  # Acceleration of the gravitational field at the stateVec location.
        for n in range(0, self.MAX_DEG+1): # Go through all the desired orders and compute the geoid corrections to the sphere.
            term = 0. # Contribution to the potential from the current degree and all corresponding orders.
            for k in range(n+1): # Go through all the orders corresponding to the currently evaluated degree.
                norm = np.sqrt((2*n+1)*m.factorial(n-k)/m.factorial(n+k))
                term += norm * as_legendre(k,n,m.cos(theta)) * (self.Cc[n][k]*m.cos( k*phi ) + self.Sc[n][k]*m.sin( k*phi ))
                    
            grav_acc += -(n+1)*m.pow(self.R/r, n+1)/r * term # Add the contribution from the current degree (derivative of the potential)
            
        grav_acc *= self.mu/self.R # Final correction.

        " Compute the acceleration due to the gravity potential at the given point. "
        grav_accel = -(r_sat/r) * abs(grav_acc)

        return grav_accel
    
    def atmosphere(self, r):
        alt = np.linalg.norm(r)-self.R
        lat, lon = self.pos2coord(r)
        return self.atm.at(self.sunclock,alt,lat,lon)

class Movody:
    '''Generic class for moving bodies = (movodies) in space'''

    def __init__(self, mass, surface, position0, velocity0, t0 = 0):
        self.m = mass
        self.S = surface
        self.r = position0.copy()
        self.v = velocity0.copy()
        self.t = t0
        self.traj = np.array([position0]) 
        self.vel = np.array([velocity0])
        self.time = [self.t]

    def clocktick(self,dt,acc):
        self.t += dt
        dv = acc*dt
        self.r += self.v*dt + dv*dt*0.5 #Euler forward integration
        self.v += dv
        self.traj = np.vstack((self.traj,self.r))
        self.vel = np.vstack((self.vel,self.v))
        self.time.append(self.t)

class Thermovody(Movody):

    def __init__(self, mass, surface, temperature, position0, velocity0, t0 = 0):
        super().__init__(mass, surface, position0, velocity0, t0)
        self.T = temperature
        self.temp = [self.T]

    def clocktick(self,dt,acc,q):
        super().clocktick(dt,acc)
        q_eff = q - cst.sigma_r*self.T**4
        self.T += q_eff*4*self.S*dt/(self.C_cal*self.m)
        if(self.T <= 0): print("Negative T : ",self.T, q_eff)
        self.temp.append(self.T)

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
        self.targets = []
        self.lostraj = {}
        self.clock = 0.
        

    def clocktick(self, dt):
        self.clock += dt
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
        v = self.orbit.getVel(self.nu)+self.orbit.rel2abs.dot(self.orbit.getVel2Rel(self.nu).dot(np.array([dvx,dvy,dvz]).T)).T
        projectile.propel(r,v,self.t)
        self.addTarget(projectile)
        Mechanics.getMech().add_animate(projectile)
        return projectile

    def addTarget(self, projectile):
        self.targets.append(projectile)
        self.tracking = True
        r = self.orbit.getPos(self.nu)
        self.linesight = projectile.r-r
        if(np.linalg.norm(self.linesight)!=0.):
            self.linesight = self.linesight/np.linalg.norm(self.linesight)
        self.lostraj[projectile] = np.array([self.linesight.copy()])
        self.traj = np.array([r.copy()])
        #self.setattitude(-self.linesight)

    def track(self):
        r = self.orbit.getPos(self.nu)
        v = self.orbit.getVel(self.nu)
        self.traj = np.vstack((self.traj,r))
        self.speed = np.vstack((self.speed,v))
        for obj in self.targets:
            self.linesight = obj.r-r
            self.linesight = self.linesight/np.linalg.norm(self.linesight)
            self.lostraj[obj] = np.vstack((self.lostraj[obj],self.linesight))
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

class Projectile(Thermovody):
    ''' Passive object representing a solid sphere and recording its kinetic and internal state trajectory '''

    smooth = 'smooth'
    coarse = 'coarse'

    def __init__(self, diameter, material, finition = coarse):
        self.d = diameter
        self.S = 0.25*m.pi*diameter**2
        self.m = m.pi/6.*diameter**3*material.mass_density
        self.surftype = finition
        self.C_cal = material.caloric_capacity
        self.k_c = material.caloric_conductivity
        self.T = 200

    def propel(self, r_ini, v_ini, t_ini = 0):
        super().__init__(self.m,self.S,self.T,r_ini,v_ini,t_ini)
        self.r_0 = r_ini
        self.r = r_ini.copy()
        self.v_0 = v_ini
        self.v = v_ini.copy()
        self.traj = np.array([r_ini])
        self.vel = np.array([v_ini])
        self.temp = [self.T]

    def clocktick(self,dt,acc,q):
        super().clocktick(dt,acc,q)
        # self.r += self.v*dt + dv*dt*0.5
        # self.v += dv
        # self.traj = np.vstack((self.traj,self.r))
        # self.vel = np.vstack((self.vel,self.v))

class Ghostectile(Thermovody):

    def __init__(self,projectile):
        self.real = projectile
        p = projectile
        super().__init__(p.m, p.S, p.T, p.r, p.v, p.t)
        self.surftype = p.surftype
        self.C_cal = p.C_cal
        self.k_c = p.k_c
        self.d = p.d
        self.r_0 = p.r.copy()
        self.v_0 = p.v.copy()



class Orbit:
    ''' Classical keplerian orbit comprising :
    * semi-major axis,
    * eccentricity,
    * inclination,
    * ascending node argument,
    * perigea argument'''

    def __init__(self,params):
        self.a, self.e, self.i, self.Om, self.w = params
        self.n = m.sqrt(cst.mu_G/self.a**3)
        self.b = self.a*m.sqrt(1-self.e**2)
        self.rel2earth = np.array([[m.cos(self.Om)*m.cos(self.w)-m.sin(self.Om)*m.cos(self.i)*m.sin(self.w), -m.sin(self.w)*m.cos(self.Om)-m.cos(self.w)*m.cos(self.i)*m.sin(self.Om), m.sin(self.i)*m.sin(self.Om)],
                        [m.cos(self.w)*m.sin(self.Om)+m.sin(self.w)*m.cos(self.i)*m.cos(self.Om), -m.sin(self.Om)*m.sin(self.w)+m.cos(self.Om)*m.cos(self.i)*m.cos(self.w), -m.sin(self.i)*m.cos(self.Om)],
                        [m.sin(self.i)*m.sin(self.w), m.sin(self.i)*m.cos(self.w), m.cos(self.i)]])
        self.earth2abs = np.array(Earth().tilt_rot.to_matrix())
        self.rel2abs = self.earth2abs.dot(self.rel2earth)


    @classmethod
    def retrieve(cls,r_ini, v_ini, r_vernal = np.array([1,0,0]), r_pole = np.array([0,0,1])):
        abs2earth = np.array((Earth().tilt_rot.inverted()).to_matrix())
        r_ini = abs2earth.dot(r_ini)
        v_ini = abs2earth.dot(v_ini)
        r_0 = np.linalg.norm(r_ini)
        v_0 = np.linalg.norm(v_ini)
        r_z = np.cross(r_ini,v_ini)
        p = np.sum(r_z**2)/cst.mu_G
        a = 1/(2/r_0 - v_0**2/cst.mu_G) # Vis-viva orbital energy equation
        e = np.sqrt(1-p/a) # Momentum to ellipticity equation
        i = abs(methu.angle_vec(np.sign(np.vdot(r_z,r_pole))*r_z,r_pole,np.cross(r_z,r_pole)))
        if(i>m.pi/2):
            i = m.pi-i
        r_node = methu.plane_intersection(r_pole,r_z)
        r_p = a*(1-e)
        if (e==0):
            theta_p = 0.
        else:
            costheta = (p/r_0-1)/e
            if(costheta>1):
                costheta = 1.
            elif(costheta<-1):
                costheta = -1.
            theta_p = np.arccos(costheta)
        rotwise = -np.sign(np.vdot(r_0,v_0))
        rotmat = np.array(Quaternion(r_z,rotwise*theta_p).to_matrix())
        r_perig = r_p * rotmat.dot(r_ini/r_0)
        Om = methu.angle_vec(r_vernal,r_node,r_pole)
        wp = methu.angle_vec(r_node,r_perig,r_pole)
        return Orbit((a,e,i,Om,wp))

    def get_nu(self,r):
        cos = (self.a*(1-self.e**2)/np.asarray(r)-1)/self.e
        ascend = (np.roll(r,1)-r >= 0)
        ascend[-1] = ascend[-2]
        cos[cos>1.] = 1.
        cos[cos<-1.] = -1.
        nu_scale = (ascend*2-1)*np.arccos(cos)
        return nu_scale

    def get_slope(self,theta):
        r = self.getPos(theta)
        v = self.getVel(theta)
        horz = v-r*np.sum(r*v,axis=1).reshape((r.shape[0],1))/np.sum(r**2,axis=1).reshape((r.shape[0],1))
        vec_ref = np.cross(horz,r)
        if(vec_ref.ndim==1):
            vec_ref = vec_ref.reshape((1,vec_ref.size))
        vec_ref = vec_ref/np.linalg.norm(vec_ref,axis=vec_ref.ndim-1).reshape((vec_ref.shape[0],1))
        slope = methu.angle_vec(horz,v,vec_ref)
        return slope

    def getPos(self, theta):
        ''' Get absolute position according to anomaly given'''
        E = self.anomaly(theta,'eccentric')
        x_plane = self.a*(np.cos(E)-self.e)
        y_plane = self.b*np.sin(E)
        z_plane = np.array([0.]*x_plane.size)
        #r = self.a*(1-self.e**2)/(1 + self.e*m.cos(theta))
        return self.rel2abs.dot(np.array([x_plane,y_plane,z_plane])).transpose()

    def getVel(self, theta):
        ''' Get absolute velocity according to anomaly and derived expressions '''
        v_x = -np.sqrt(cst.mu_G/(self.a*(1-self.e**2)))*np.sin(theta)
        v_y = np.sqrt(cst.mu_G/(self.a*(1-self.e**2)))*(self.e+np.cos(theta))
        v_z = np.array([0.]*v_x.size)
        return self.rel2abs.dot(np.array([v_x,v_y,v_z])).transpose()

    def getVel2(self,theta):
        ''' Get absolute velocity according to energy relation and tangent vector'''
        r = self.a*(1-self.e**2)/(1+self.e*m.cos(theta))
        vnorm = m.sqrt(cst.mu_G*(2/r-1/self.a))
        v = vnorm*np.array([1,0,0])
        return self.getVel2Rel(theta).dot(v.T)

    def getVel2Rel(self, theta):
        ''' Get transformation matrix from the Fresnel frame to the Relative frame'''
        t_vec = np.array([-self.a*m.sin(theta),self.b*m.cos(theta),0])
        n_vec = np.array([-self.b*m.cos(theta),-self.a*m.sin(theta),0])
        t_vec = t_vec/np.linalg.norm(t_vec)
        n_vec = n_vec/np.linalg.norm(n_vec)
        return np.array([t_vec,n_vec,np.array([0,0,1])]).transpose()

    def getRad2Rel(self,theta):
        ''' Get transformation matrix from a the Polar frame to the Relative frame'''
        r_vec = np.array([m.cos(theta),m.sin(theta),0])
        pol_vec = np.array([-m.sin(theta),m.cos(theta),0])
        return np.array([r_vec,pol_vec,np.array([0,0,1])]).transpose()

    def anomaly(self, anom, output = 'true'):
        ''' Convert an anomaly type into a different one'''
        res_anomaly = anom
        corr = anom%(2*m.pi)>m.pi
        if(output=='true'):
            f = lambda x : x - self.e*np.sin(x) - anom
            ec_anomaly = sp.optimize.root(f,anom,method='lm').x[0]
            res_anomaly = (1-2*corr)*np.arccos( (np.cos(ec_anomaly)-self.e) / (1 - self.e*np.cos(ec_anomaly))) + 2*m.pi*corr
        elif(output=='mean'):
            ec_anomaly = (1-2*corr)*np.arccos( (self.e+np.cos(anom)) / (1+self.e*np.cos(anom))) + 2*m.pi*corr # eccentric anomaly
            res_anomaly = ec_anomaly - self.e*np.sin(ec_anomaly)
        elif(output=='eccentric'):
            ec_anomaly = (1-2*corr)*np.arccos( (self.e+np.cos(anom)) / (1+self.e*np.cos(anom))) + 2*m.pi*corr # eccentric anomaly
            res_anomaly = ec_anomaly
        return res_anomaly

class Atmosphere:
    ''' Based on the NASA MSISE-90 model implementation, yields a thermodynamic stated air according to the specified time and location'''

    def __init__(self, year = 2000, doy = 1, model = 'advanced'):
        self.model = model
        self.year = year
        self.doy = doy
        self.rho_hist = []
        self.flags = atm.nrlmsise_flags()
        self.flags.switches[:] = [1]*(len(self.flags.switches)) # output in kg/m³

    def at(self, timesec, altitude, latitude, longitude):
        ''' Air properties at a given position
        ----
        Params : 
        * timesec : time of the day (seconds)
        * altitude : (meters)
        * latitude : (radians)
        * longitude : (radians)'''
        if(self.model=='advanced'):
            output = atm.nrlmsise_output()
            inputatm = atm.nrlmsise_input(year = self.year, doy = self.doy, sec = timesec, alt=altitude*0.001, g_lat=latitude, g_long=longitude, lst=0.0, f107A=150., f107=150.0, ap=4.0, ap_a=None)
            atm.lstCalc(inputatm)
            gtd7(inputatm,self.flags,output)
            rho = output.d[5]
            T = output.t[1]
        else:
            rho = density(altitude*0.001)
            T = 214.
        air = Atmosphere.Air(rho,T)
        # self.rho_hist.append(rho)
        return air

    def profile(self, hscale=range(0,400,1), latitude = 0., longitude = 0.):
        rho = []
        o = atm.nrlmsise_output()
        i = atm.nrlmsise_input(year = self.year, doy = self.doy, sec = 0., g_lat=latitude, g_long=longitude, lst=0.0, f107A=150., f107=150.0, ap=4.0, ap_a=None)
        for h in hscale:
            i.alt = h
            gtd7(i,self.flags,o)
            rho.append(o.d[5])
        return np.array(rho), hscale

    class Air:

        N_A = cst.N_A # Avogadro
        MM_air = cst.MM_air # kg/mol
        R_m = cst.R/MM_air # Jmol/kgK
        gamma = cst.gamma # gas caloric constant
        sigma_c = cst.sigma_c # m² mean collisionnal cross section
        b_air = cst.b_air # kg/m.s(K)^0.5
        S_air = cst.S_air # K

        def __init__(self,rho,T):
            self.rho = rho # mass density
            self.nv = rho*self.N_A/self.MM_air # particle density
            self.mean_free_path = 1/(m.sqrt(2)*self.nv*self.sigma_c) # mean free path between collisions
            self.T = T # kinetic air temperature
            self.p = self.rho*self.R_m*self.T
            self.c_a = m.sqrt(self.gamma*self.R_m*self.T)
            self.visc = self.b_air*self.T**(1.5)/(self.T+self.S_air) # Sutherland law for gases : mass viscosity
            


class DynamicBox:
    ''' Helper class based on the Decorator Design Pattern for implementing dynamics on objects'''

    def __init__(self, earth, time = 0):
        self.earth = earth
        self.atmosphere = earth.atm
        self.timeline = [time]
        self.clock = time

    def set_movody(self,movody):
        self.movody = movody

    def clocktick(self,dt):
        self.clock += dt
        self.timeline.append(self.clock)


    def grav(self):
        a_p = self.earth.gravity(self.movody.r)
        # a_p2 = -self.earth.mu * self.movody.r/np.linalg.norm(self.movody.r)**3
        # print("Precise : ",a_p1, "\n", "Coarse : ", a_p2)
        return a_p

    def N_coll(self,n_mean):
        n_coll = n_mean
        if(n_mean <= 1000):
            n_coll = int(np.round((np.random.normal(n_mean,1/m.sqrt(n_mean)))))
        return n_coll

    def drag(self,dt):
        air = self.earth.atmosphere(self.movody.r)

        volume = self.movody.S*np.linalg.norm(self.movody.v)*dt # Tube volume covered by projectile during dt
        

        n_volmean = air.nv
        # n_volmean = density((np.linalg.norm(obj.r)-self.earth.R)*0.001)*N_A/MM_air
        lpm = 1/(m.sqrt(2)*n_volmean*air.sigma_c) # mean free path formulae
        n_mean = n_volmean*volume
        n_coll = self.N_coll(n_mean)
        dv_coll = - n_coll*air.MM_air/(air.N_A*self.movody.m) * self.movody.v
        
        # a_p = -self.movody.r/np.linalg.norm(self.movody.r)*cst.mu_G/(np.linalg.norm(self.movody.r))**2

        dv_p = dv_coll
        
        return dv_p, air

class SpaceBox(DynamicBox):
    ''' Helper class that supervises the dynamics of the satellite during its orbiting around Earth'''

    def __init__(self, earth, time = 0):
        super().__init__(earth,time)

    def clocktick(self, dt):
        super().clocktick(dt)
        dv = self.drag(dt)
        self.sat.clocktick(dt)


    def set_satellite(self,satellite):
        self.sat = satellite

class GhostsonicBox(DynamicBox):

    def __init__(self, earth, time = 0):
        super().__init__(earth,time)
        self.sigma_d = 0
        self.pathmark = []
        self.machline = []
        self.rholine = []
        self.presline = []
        self.Knline = []
        self.Reline = []
        self.dragline = []
        self.heatline = []
        self.markindex = []

    def clocktick(self,dt):
        super().clocktick(dt)
        #dv_d, air_d = super().drag(dt)
        a_d, air, C_D, q, Kn_mix = self.drag()
        a_g = self.grav()
        self.movody.clocktick(dt,a_d+a_g,q)
        self.rholine.append(air.rho)
        self.Knline.append(Kn_mix)
        v_p = np.linalg.norm(self.movody.v)
        self.machline.append(v_p/air.c_a)
        self.presline.append(air.p)
        self.Reline.append(air.rho*v_p*self.movody.d/air.visc)
        self.dragline.append(C_D)
        self.heatline.append(q)
        if self.differs(): 
            self.record(air)
            self.request_sim()

    def differs(self):
        alt_now = (np.linalg.norm(self.movody.r)-self.earth.R)
        alt_prev = np.linalg.norm(self.pathmark[-1][0])-self.earth.R
        # print(alt_prev)
        # print("Alt Difference : ",(alt_now-alt_prev)/alt_prev)
        return abs((alt_now-alt_prev)/alt_prev) >= 0.1 # Altitude change of at least 10 %


    def record(self,air):
        self.markindex.append(max(len(self.rholine)-1,0))
        self.pathmark.append((self.movody.r.copy(),self.movody.v.copy(),self.movody.T,air))

    def request_sim(self):
        self.dsmc.add_sim(*self.pathmark[-1])

    def set_movody(self,movody):
        super().set_movody(movody)
        movody.aeroreg = AeroRegime.Free_Molecular_Flow
        # self.sigma_d = HypersonicBox._surftype[movody.surftype]
        self.regime = FMF(movody)
        self.record(self.earth.atmosphere(self.movody.r))

    def set(self,dsmc):
        self.dsmc = dsmc
        self.request_sim()

    def drag(self):
        air = self.earth.atmosphere(self.movody.r)
        v_p = np.linalg.norm(self.movody.v)
        Kn_mix = self.find_regime()
        C_D = self.regime.get_dragC(air)
        q = self.regime.get_heat(air)
        f_drag = - 0.5*C_D*self.movody.S*air.rho*v_p*self.movody.v/self.movody.m
        if(q<=0): print("Negative Heat in drag : ", q)
        return f_drag, air, C_D, q, Kn_mix

    def find_regime(self):
        air = self.earth.atmosphere(self.movody.r)
        v_p = np.linalg.norm(self.movody.v)
        # if(self.movody.T<0): print("Negative T ! : ",self.movody.temp)
        v_therm = np.sqrt(8*air.R_m*self.movody.T/np.sqrt(m.pi))
        mfp_mix = v_therm/(air.nv*v_p*air.sigma_c) # lambda_m = v_moy/(n*U*sigma)
        Kn_mix = mfp_mix/self.movody.d
        return Kn_mix

# Strategy Desgin pattern for dividing the different flow regimes in the box.
class HypersonicBox(DynamicBox):
    ''' Helper class that supervises the dynamics of the projectile during its descent in the atmosphere'''

    def __init__(self, earth, time = 0):
        super().__init__(earth,time)
        # self.sigma_d = 0
        self.machline = []
        self.rholine = []
        self.presline = []
        self.Knline = []
        self.Reline = []
        self.dragline = []
        self.heatline = []
        self.regimes = {}
        self.ghostbox = None
        self.ghost = None

    def clocktick(self,dt):
        super().clocktick(dt)
        #dv_d, air_d = super().drag(dt)
        a_d, air, C_D, q, Kn_mix = self.drag()
        a_g = self.grav()
        self.movody.clocktick(dt,a_d+a_g,q)
        self.rholine.append(air.rho)
        self.Knline.append(Kn_mix)
        v_p = np.linalg.norm(self.movody.v)
        self.machline.append(v_p/air.c_a)
        self.presline.append(air.p)
        self.Reline.append(air.rho*v_p*self.movody.d/air.visc)
        self.dragline.append(C_D)
        self.heatline.append(q)

    def set_movody(self,movody):
        super().set_movody(movody)
        movody.aeroreg = AeroRegime.Free_Molecular_Flow
        self.dsmc = DSMC(movody)
        # self.sigma_d = HypersonicBox._surftype[movody.surftype]
        self.regimes[AeroRegime.Free_Molecular_Flow] = FMF(movody)
        self.regimes[AeroRegime.Transitional_Flow] = TF(movody,self.dsmc)


    def dragsim(self,dt):
        volume = self.movody.S*np.linalg.norm(self.movody.v)*dt # Tube volume covered by projectile during dt

        air = self.earth.atmosphere(self.movody.r)
        n_volmean = air.nv
        n_mean = n_volmean*volume
        n_coll = self.N_coll(n_mean)
        
        dv_coll = - n_coll*air.MM_air/(air.N_A*self.movody.m) * self.movody.v
    
        a_p = -self.movody.r/np.linalg.norm(self.movody.r)*cst.cst.mu_G/(np.linalg.norm(self.movody.r))**2

        dv_p = a_p*dt + dv_coll
        return dv_p, air

    def drag(self):
        air = self.earth.atmosphere(self.movody.r)
        v_p = np.linalg.norm(self.movody.v)
        Kn_mix = self.find_regime()
        C_D = self.regime.get_dragC(air)
        q = self.regime.get_heat(air)
        f_drag = - 0.5*C_D*self.movody.S*air.rho*v_p*self.movody.v/self.movody.m

        return f_drag, air, C_D, q, Kn_mix

    def find_regime(self):
        air = self.earth.atmosphere(self.movody.r)
        v_p = np.linalg.norm(self.movody.v)
        v_therm = np.sqrt(8*air.R_m*self.movody.T/np.sqrt(m.pi))
        mfp_mix = v_therm/(air.nv*v_p*air.sigma_c) # lambda_m = v_moy/(n*U*sigma)
        Kn_mix = mfp_mix/self.movody.d
        if Kn_mix <= 15 and (self.ghostbox is None):
            self.ghost = Ghostectile(self.movody)
            self.namecode = 'Ghost Projectile Cloned'
            Mechanics().schedule.urge(self.namecode, self.ghostready, delay = 0., duration = TimeGenerators().ghostectileTimeGenerator(self.ghost),severed = True)
        if(Kn_mix >= 10):
            self.regime = self.regimes[AeroRegime.Free_Molecular_Flow]
            self.movody.aeroreg = AeroRegime.Free_Molecular_Flow
        else:
            self.regime = self.regimes[AeroRegime.Transitional_Flow]
            self.movody.aeroreg = AeroRegime.Transitional_Flow
        return Kn_mix

    def ghostready(self):
        self.ghostbox = Mechanics().add_animate(self.ghost,self.namecode)
        self.ghostbox.set(self.dsmc)

class AeroRegime:

    Free_Molecular_Flow = 'fmf'
    Transitional_Flow = 'tf'

    _surftype = {Projectile.smooth: 0, Projectile.coarse : 1}

    def __init__(self, thermovody):
        self.thermovody = thermovody
        self.sigma_d = AeroRegime._surftype[thermovody.surftype]

    def get_dragC(self,air):
        return 2

    def get_heat(self,air):
        pass



class FMF(AeroRegime):

    def get_dragC(self,air):
        v_p = np.linalg.norm(self.thermovody.v)
        S = np.sqrt(air.gamma/2)*v_p/air.c_a
        C_D = m.e**(-S**2)/(np.sqrt(2)*S**3)*(1+2*S**2) + (4*S**4+4*S**2-1)/(2*S**4)*m.erf(S) + 2*self.sigma_d*np.sqrt(np.pi)/(3*S)
        return C_D

    def get_heat(self,air):
        v_p = np.linalg.norm(self.thermovody.v)
        S = np.sqrt(air.gamma/2)*v_p/air.c_a
        e_moy = cst.k_B * air.T * 0.5 # kT/2
        N_i = air.nv * np.sqrt(air.R_m*air.T/(2*m.pi)) # n* sqrt(kT/2pi m)
        diff = self.sigma_d
        Tr = self.thermovody.T
        zeta = 5
        q = diff*N_i*e_moy* ((2*S**2+1+(4+zeta)*(1-Tr/air.T))*m.sqrt(m.pi)*(m.erf(S)/(2*S)*(1.5+S)+m.e**(-S**2)/(2*m.sqrt(m.pi)))-m.erf(S)/S)
        # print("FMF Heat : ", q, " | FMF proj T : ", air.T)
        return q

class TF(AeroRegime):

    def __init__(self,thermovody,dsmc):
        super().__init__(thermovody)
        self.dsmc = dsmc
        self.dragCall = False
        self.heatCall = False
#        self.C_D
#        self.q

    def get_dragC(self,air):
        self.dragCall = not self.dragCall
        if(self.dragCall != self.heatCall):
            self.C_D, self.q = self.dsmc.exertion(air)
        return self.C_D

    def get_heat(self,air):
        self.heatCall = not self.heatCall
        if(self.dragCall != self.heatCall):
            self.C_D, self.q = self.dsmc.exertion(air)
        return self.q

class TimeGenerators(metaclass=Singleton):

    def __init__(self):
        self.earth = Earth()
        self.atm = self.earth.atm

    def get_dt(self,projectile):
        if(projectile.aeroreg=='fmf'):
            dt = 10**(-3)*self.earth.R/np.linalg.norm(projectile.v)
        else:
            air = self.earth.atmosphere(projectile.r)
            v_p = np.linalg.norm(projectile.v)
            t_coll = 1/(air.nv*v_p*air.sigma_c)
            dt = t_coll/3
        return dt

    def ghostectileTimeGenerator(self,ghostectile, N_itmax = 10**4):
        alt = (np.linalg.norm(ghostectile.r)-self.earth.R)*0.001
        n_it = 0
        while (alt >= 50 and n_it < N_itmax and alt>=0):
            alt = (np.linalg.norm(ghostectile.r)-self.earth.R)*0.001
            dt = 10**(-3)*self.earth.R/np.linalg.norm(ghostectile.v)
            n_it += 1
            yield dt

    def projectileTimeGenerator(self,projectile, N_itmax = 10**4):
        alt = (np.linalg.norm(projectile.r)-self.earth.R)*0.001
        n_it = 0
        while (alt >= 50 and n_it < N_itmax and alt>=0):
            alt = (np.linalg.norm(projectile.r)-self.earth.R)*0.001
            dt = self.get_dt(projectile)
            # print("yielding proj: ",dt)
            n_it += 1
            yield dt

    def orbitalTimeGenerator(self,movody, duration, N_itmax = 10**4):
        dt = 10**(-3)*self.earth.R/np.lin*alg.norm(projectile.v_0)
        t = 0
        while(t<duration and n_it < N_itmax):
            t += dt
            if(t>duration):
                dt -= t-duration
            yield dt
 
class Mechanics(metaclass=Singleton):
    ''' Master class for executing events, ticking the clock and notifying animated objects'''
    boxing = {Projectile: HypersonicBox, Ghostectile: GhostsonicBox}

    def __init__(self):
        self.animate = []
        print('new Mech :', self)
        self.clock = 0.
        self.mainline = Mechanics.Timeline(self.clock,self)
        self.timelines = [self.mainline]
        self.linekeys = {0:0}
        self.timegens = []
        self._timeline = [self.clock]
        self.boxes = {}

    @classmethod
    def initialise(cls):
        ''' Method for replacing the old mechanics object by a new one, must by imperatively called first if program in a script'''
        cls._instances.pop(cls,'empty')
        return Mechanics()

    @classmethod
    def getMech(cls):
        ''' Convenient method for accessing the single mechanics instance'''
        return Mechanics._instances[cls]

    @property
    def timeline(self):
        return self.mainline.timeline
        # return self._timeline


    def set(self,earth,schedule):
        self.earth = earth
        self.schedule = schedule

    def add_animate(self, obj, key = 0):
        if(obj.__class__ in Mechanics.boxing):
            box = Mechanics.boxing[obj.__class__](self.earth,self.clock)
            box.set_movody(obj)
            self.boxes[obj] = box   
            obj = box
            # self.hypersonicbox = HypersonicBox(self.earth,self.clock)
            # self.hypersonicbox.set_airobj(obj)
            # obj = self.hypersonicbox           
        tl = self.timelines[self.linekeys[key]]
        tl.animate.append(obj)
        self.animate.append(obj)
        return obj

    def fetch_events(self):
        evts = self.schedule.at(self.clock)
        for evt in evts:
            # print('Fetched : '+ evt.name)
            tgen = evt.timegen
            if(evt.severed): 
                self.linekeys[evt.name] = len(self.timelines)
                newtl = Mechanics.Timeline(self.clock)
                newtl.timegens.append(tgen)
                self.timelines.append(newtl) 
            evt.trigger()
            # print('Generator : ',tgen)
            if(evt.severed):               
                newtl.start()
            else:
                self.mainline.timegens.append(tgen)
                self.timegens.append(tgen)

    def notify(self, dt):
        self.clock += dt
        self.fetch_events()
        return 'ok'

    def alert(self):
        time = self.schedule.next()
        if time: 
            self.clock = time
            self.fetch_events()
        return (time != None)          

    def start(self):
        # 1 - list of events on schedule
        # 2 - calculation of time pace
        # 3 - update of animate objects
        self.fetch_events()
        # Loop until job finished
        self.mainline.run()
        
        # while self.mainline.is_alive():
        #     pass #sleep
        # while len(self.timegens)!=0:
        #     dts = []
        #     for i,tgen in enumerate(self.timegens.copy()):
        #         try:
        #             dts.append(next(tgen))
        #         except StopIteration:
        #             self.timegens.remove(tgen)
        #             print('depleted')
        #         finally:
        #             if(len(self.timegens)==0):
        #                 print('Mechanics stopped')
        #                 return
        #     dt = min(dts)
        #     self.move(dt)
        #     self.clock += dt
        #     self.timeline.append(self.clock)
        #     self.fetch_events()

    # def move(self, dt):
    #     if(dt!=0):
    #         for obj in self.animate:
    #                 obj.clocktick(dt)

    class Timeline(Thread):

        def __init__(self,clockset = 0, call = None):
            super().__init__()
            self.clock = clockset
            self.timeline = [self.clock]
            self.animate = []
            self.timegens = []
            self.history = {}
            self.call = call

        def run(self):
#            super().run()
            # Loop until job finished
            while len(self.timegens)!=0:
                loop = False
                dts = []
                for i,tgen in enumerate(self.timegens.copy()):
                    try:
                        dts.append(next(tgen))
                    except StopIteration:
                        self.timegens.remove(tgen)
                        print('depleted')
                        if(len(self.timegens)==0):
                            loop = False
                            if self.call:
                                print('Calling back, last chance')
                                loop = self.call.alert()
                            if not loop :
                                print('No more events')
                                return
                if loop: print('looping'); continue
                dt = min(dts)
                self.move(dt)
                self.clock += dt
                self.timeline.append(self.clock)
                if self.call :
                    self.call.notify(dt)

        def move(self, dt):
            if(dt!=0):
                for obj in self.animate:
                        obj.clocktick(dt)

    class SynClock(float):

        def __init__(self,number = 0.):
            self._value = number

        def __iadd__(self,other):
            self._value += other

        @property
        def value(self):
            return self._value


class Schedule:
    ''' Basic implementation of a schedule containing events disposed on a timetable. 
    Schedule can be consulted at any time for fetching events on schedule for one time only read.'''
    
    def __init__(self):
        self.events = {}
        self.times = []
    
    def plan(self, eventname, func, *args, time = 0, duration = iter([0]), severed = False):
        evt = Event(eventname,duration,func,*args, severed = severed)
        if(time in self.events):
            self.events[time].append(evt)
        else:
            self.events[time] = [evt]
            self.times.append(time)
        self.times.sort()

    def urge(self, eventname, func, *args, delay = 0, duration = iter([0]), severed = False):
        self.plan(eventname, func, *args, time = self.record+delay, duration = duration, severed = severed)

    def at(self, time):
        self.record = time
        if(len(self.times)!=0 and time>=self.times[0]):
            return self.events[self.times.pop(0)]
        else:
            return []

    def next(self):
        if len(self.times)!=0: return self.times[0]
        return None

class Event:
    ''' Object containing a name, a task function and a slicing time generator'''

    def __init__(self, name, timegen, func, *args, **kwargs):
        self.name = name
        self.func = func
        self.args = args
        self.timegen = timegen
        self.severed = kwargs.pop("severed",False)

    def trigger(self):
        print(self.name)
        if(self.func is not None):
            self.func(*self.args)   

class Reference:
    ''' Convenient class for regrouping time and frame references of the code '''
    _refdate = dt.datetime(2000,1,1,12)

    @classmethod
    def get_refdate(cls):
        return Reference._refdate

    @classmethod
    def get_refloc(cls):
        return 35.6895*m.pi/180, 139.69171*m.pi/180, 'Tokyo' # Tokyo coordinates : Latitude, Longitude

class DSMC:

    sparta_folder = "dsmc"
    batch_file_name = "dsmc_cmd.bat"
    cygwin_path = "C:\\cygwin\\bin\\mintty.exe"
    bash_path = "C:\\cygwin\\bin\\bash.exe"
    sparta_instr = "spa_serial < "
    sparta_exe = os.getcwd()+"\\"+sparta_folder+"\\"+"spa_serial"
    MAX_SIM_PARA = 4

    def __init__(self, thermovody):
        self.thermovody = thermovody
        self.path2dsmc = os.getcwd()+"\\"+self.sparta_folder+"\\"
        filename = self.path2dsmc + self.batch_file_name
        instr = self.prepare_batch_file(filename)
        self.regimes = {}
        self.simlist = []
        self.running = {}
        self.finished = {}
        self.simres = {}
        self.simind = 0
        self.writer = dsmc.Sparta_Writer()
        self.sentinel = Thread(target = self.watch)
        # self.reader = dsmc.Sparta_Reader()
        # self.launched = False
        # self.launch(filename,instr)

    def watch(self):
        self.watching = True
        self.delay = 5 # s
        while self.watching:
            for fileid,proc in self.running.copy().items():
                if proc.poll():
                    print(fileid," Simulation results : ",proc.communicate())
                    self.running.pop(fileid)
                    self.retrieve_results(fileid)
                    self.finished[fileid] = proc
                    if(len(self.finished.values()) == len(self.simlist)) : self.watching = False
            time.sleep(self.delay)      

    def prepare_batch_file(self,filename):
        f = pathlib.Path(filename)
        instr = "spa_serial < in.axinit"
        # if(not f.exists()):
        cd = "cd "+ "\""+str(pathlib.Path(self.path2dsmc).absolute())+"\""
        cmd = "@"+self.bash_path + " --login -c " + cd + " " + instr
        f.write_text(cmd)
        return shlex.split(instr)

    def launch(self,filename,instr):
        if(not self.sentinel.is_alive()): self.sentinel.start()
        print("Launching Simulation : "+filename.split("\\")[-1])
        sim = subprocess.Popen(instr,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = self.sparta_folder, universal_newlines = True, shell=False)
        self.running[filename] = sim
        self.simind += 1
        return sim
        # self.proc.wait()
        # print(self.proc.communicate())

    def add_sim(self,r,v,T,air):
        alt = (np.linalg.norm(r)-Earth.R)*0.001
        speed = np.linalg.norm(v)
        sim_params = self.prepare_sim(alt,speed,T,air)
        fresh = len(self.simlist)==0
        simname = self.writer.write_simulation(sim_params,fresh)
        self.simlist.append(sim_params)
        print("Preparing Simulation : ",simname.split("\\")[-1])
        self.regimes[simname] = (DSMC.Regime(r,v,T,air,self.thermovody))
        if(fresh):
            self.launch(simname,[self.sparta_exe,"-in",simname])
        elif(len(self.running.values()) < self.MAX_SIM_PARA and len(self.finished.values()) == 1):
            # self.launch(simname,self.sparta_instr+simname)
            self.launch(simname,[self.sparta_exe,"-in",simname])



    def prepare_sim(self,alt,speed,T_proj,air):
        Ma = speed/air.c_a
        param = {"name":"{:.2f}km_{:d}".format(alt,int(Ma))}
        param[dsmc.SimInput.AIR_NDENSITY] = air.nv
        param[dsmc.SimInput.AIR_TEMPERATURE] = air.T
        param[dsmc.SimInput.SPEED] = speed
        param[dsmc.SimInput.PROJ_TEMPERATURE] = T_proj
        return param

    def retrieve_results(self,file):
        # output = self.reader.read(file)
        self.regimes[file].C_D = 2+0.1*(len(self.finished.values()))
        self.regimes[file].q = 0.
        print("Retrieving from ",file)
        while(len(self.running.values()) < self.MAX_SIM_PARA):
            simname = self.simlist[self.simind]["simname"]
            self.launch(simname,[self.sparta_exe,"-in",simname])
        pass

    def exertion(self,air):
        reg = DSMC.Regime(self.thermovody.r,self.thermovody.v,self.thermovody.T,air,self.thermovody)
        keyids = self.find_closest(reg)
        while not all(id in self.finished for id in keyids):
            print("Waiting for simulation at {:d} km".format(int((np.linalg.norm(self.thermovody.r)-Earth.R)*0.001)))
            time.sleep(3) # waiting for simulation
        return self.combine(reg,[self.regimes[id] for id in keyids])

    def combine(self,reg,basis):
        dist = np.array([])
        drags = []
        heats = []
        for r in basis:
            dist = np.hstack((dist,reg.distance(r)))
            drags.append(r.C_D)
            heats.append(r.q)
        tot = sum(dist)
        dist = dist/tot
        C_D = np.array(drags).dot(dist)
        heat = np.array(heats).dot(dist)
        return C_D, heat


    def find_closest(self,reg):
        dist = {}
        for id,r in self.regimes.items():
            d = reg.distance(r)
            dist[id] = d
        try:
            ids = methu.minima(dist,2)
        except(Exception):
            ids = self.finished.keys()
        return ids

    def abort_all(self):
        for f,p in self.running.items():
            p.terminate()
        self.watching = False


    class Regime:

        Kn_w = 0.5
        Ma_w = 0.5

        def __init__(self,r,v,T,air,thermovody):
            self.air = air
            self.U = np.linalg.norm(v)
            self.L = thermovody.d
            self.alt = np.linalg.norm(r)-Earth.R
            self.Ma = self.U/self.air.c_a
            v_therm = np.sqrt(8*air.R_m*T/np.sqrt(m.pi))
            self.Kn = self.air.mean_free_path/self.L*v_therm*np.sqrt(2)/self.U

        def distance(self,reg):
            return np.sqrt((self.Ma-reg.Ma)**2*self.Ma_w + (self.Kn-reg.Kn)**2*self.Kn_w)
            
