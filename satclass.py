# -*- coding: utf-8 -*-
"""
Created on Thu May 25 11:42:35 2017

@author: ASiapan
"""

import math as m
import numpy as np
import scipy as sp
import csv
import atmosphere.nrlmsise_00_header as atm
from atmosphere.nrlmsise_00 import gtd7
import datetime as dt
from density import density
from mathutils import Quaternion
#from ALE import mu_G, R_G

global mu_G, R_G
mu_G = 398600 * 10**9 # m³/s²
R_G = 6378 * 10**3 # m 


MM_air = 0.0289644 # kg/mol
N_A = 6.023 * 10**23
sigma_c = 1e-19 # m²

class Satellite:
    
    def __init__(self, orbit, nu_ini, width):
        self.orbit = orbit
        self.tracking = False
        self.nu_0 = self.nu = nu_ini
        self.M_0 = self.M = orbit.anomaly(nu_ini,'mean')
        self.traj = np.array([orbit.getPos(nu_ini)])
        self.speed = np.array([orbit.getVel(nu_ini)])
        self.width = width
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
        v = self.orbit.getVel(self.nu)+self.orbit.rel2abs.dot(self.orbit.getVel2Rel(self.nu).dot(np.array([dvx,dvy,dvz]).T)).T
        projectile.propel(r,v)
        self.setTarget(projectile)
        Mechanics().add_animate(projectile)
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

class Orbit:

    def __init__(self,params):
        self.a, self.e, self.i, self.Om, self.w = params
        self.n = m.sqrt(mu_G/self.a**3)
        self.b = self.a*m.sqrt(1-self.e**2)
        self.rel2abs = np.array([[m.cos(self.Om)*m.cos(self.w)-m.sin(self.Om)*m.cos(self.i)*m.sin(self.w), -m.sin(self.w)*m.cos(self.Om)-m.cos(self.w)*m.cos(self.i)*m.sin(self.Om), m.sin(self.i)*m.sin(self.Om)],
                        [m.cos(self.w)*m.sin(self.Om)+m.sin(self.w)*m.cos(self.i)*m.cos(self.Om), -m.sin(self.Om)*m.sin(self.w)+m.cos(self.Om)*m.cos(self.i)*m.cos(self.w), -m.sin(self.i)*m.cos(self.Om)],
                        [m.sin(self.i)*m.sin(self.w), m.sin(self.i)*m.cos(self.w), m.cos(self.i)]])

    def getPos(self, theta):
        E = self.anomaly(theta,'eccentric')
        x_plane = self.a*(m.cos(E)-self.e)
        y_plane = self.b*m.sin(E)
        #r = self.a*(1-self.e**2)/(1 + self.e*m.cos(theta))
        return self.rel2abs.dot(np.array([x_plane,y_plane,0.]).transpose())

    def getVel(self, theta):
        v_x = -m.sqrt(mu_G/(self.a*(1-self.e**2)))*m.sin(theta)
        v_y = m.sqrt(mu_G/(self.a*(1-self.e**2)))*(self.e+m.cos(theta))
        return self.rel2abs.dot(np.array([v_x,v_y,0.])).transpose()

    def getVel2(self,theta):
        r = self.a*(1-self.e**2)/(1+self.e*m.cos(theta))
        vnorm = m.sqrt(mu_G*(2/r-1/self.a))
        v = vnorm*np.array([1,0,0])
        return self.getVel2Rel(theta).dot(v.T)

    def getVel2Rel(self, theta):
        t_vec = np.array([-self.a*m.sin(theta),self.b*m.cos(theta),0])
        n_vec = np.array([-self.b*m.cos(theta),-self.a*m.sin(theta),0])
        t_vec = t_vec/np.linalg.norm(t_vec)
        n_vec = n_vec/np.linalg.norm(n_vec)
        return np.array([t_vec,n_vec,np.array([0,0,1])]).transpose()

    def getRad2Rel(self,theta):
        r_vec = np.array([m.cos(theta),m.sin(theta),0])
        pol_vec = np.array([-m.sin(theta),m.cos(theta),0])
        return np.array([r_vec,pol_vec,np.array([0,0,1])]).transpose()

    def anomaly(self, anom, output = 'true'):
        res_anomaly = anom
        corr = (anom>m.pi)
        if(output=='true'):
            f = lambda x : x - self.e*m.sin(x) - anom
            ec_anomaly = sp.optimize.root(f,anom,method='lm').x[0]
            res_anomaly = (1-2*corr)*m.acos( (m.cos(ec_anomaly)-self.e) / (1 - self.e*m.cos(ec_anomaly))) + 2*m.pi*corr
        elif(output=='mean'):
            ec_anomaly = (1-2*corr)*m.acos( (self.e+m.cos(anom)) / (1+self.e*m.cos(anom))) + 2*m.pi*corr # eccentric anomaly
            res_anomaly = ec_anomaly - self.e*m.sin(ec_anomaly)
        elif(output=='eccentric'):
            ec_anomaly = (1-2*corr)*m.acos( (self.e+m.cos(anom)) / (1+self.e*m.cos(anom))) + 2*m.pi*corr # eccentric anomaly
            res_anomaly = ec_anomaly
        return res_anomaly

class Earth:

    def __init__(self):
        self.R = R_G
        self.mu = mu_G
        self.sidereal = 86164.1 # sidereal day duration
        self.solarday = 86400 # solar day duration
        self.w_rot = 2*m.pi/self.sidereal
        self.clock = 0.
        self.sunclock = 0.
        self.rotang = 0.
        self.phi = [self.rotang]
        self.Cc, self.Sc = self.unpackgravitymodel()
        self.MAX_DEG = 10

    def setDate(self,date):
        self.clock = (date - Reference.get_refdate()).total_seconds()%self.sidereal
        self.sunclock = (date - Reference.get_refdate()).total_seconds()%self.solarday
        self.rotang = self.clock*self.w_rot
        self.phi = [self.rotang]
        self.atm = Atmosphere(date.year, date.day, model = 'advanced')

    def clocktick(self,dt):
        self.clock += dt
        self.sunclock += dt
        self.rotang = (self.clock%self.sidereal)*self.w_rot
        self.phi.append(self.rotang)

    def pos2coord(self,r):
        g_rot = np.array(Quaternion([0,0,1],-self.rotang).to_matrix())
        x_rel, y_rel, z_rel = g_rot.dot(r)
        rnorm = np.linalg.norm(r)

        lon = m.atan2(y_rel,x_rel)
        lat = m.atan2(z_rel,m.sqrt(y_rel**2+x_rel**2))

        return lat*180/m.pi, lon*180/m.pi

    def unpackgravitymodel(self):
        """ Read the EGM96 gravitational field model coefficients from EGM96coefficients
        file and parse them to be used with computeGravitationalPotential functions.
        
        Returns
        ----------
        2-tuple of the C and S coefficients of EGM96 model. They are stored in dictionaries
            of list. The keys are degrees of the potential expansion and the values
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
                    Ccoeffs[newdeg].append([float(row[3])])
                    Scoeffs[newdeg].append([float(row[4])])
                olddeg = newdeg

        return Ccoeffs, Scoeffs

    def gravity(self,r_sat,epoch):
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
        r = np.linalg.norm(r_sat)
        theta = m.acos(r_sat[2]/r)
        phi = m.atan(r_sat[1]/r_sat[0])
        as_legendre = sp.special.lpmv
        " Find the gravitational potential at the desired point. "
        grav_pot = 0.0 # Potential of the gravitational field at the stateVec location.
        for n in range(0, self.MAX_DEG+1): # Go through all the desired orders and compute the geoid corrections to the sphere.
            term = 0. # Contribution to the potential from the current degree and all corresponding orders.
            for k in range(n+1): # Go through all the orders corresponding to the currently evaluated degree.
                term += as_legendre(k,n,m.cos(theta)) * (self.Cc[n][k]*m.cos(k*phi) + self.Sc[n][k]*m.sin( k*phi ))
                    
            grav_pot += m.pow(self.R/r, n) * term # Add the contribution from the current degree.
            
        grav_pot *= self.mu/self.R # Final correction.

        " Compute the acceleration due to the gravity potential at the given point. "
        grav_accel = -(r_sat/r) * grav_pot/r

        return grav_accel
    
    def atmosphere(self, r):
        alt = np.linalg.norm(r)-self.R
        lat, lon = self.pos2coord(r)
        return self.atm.at(self.sunclock,alt,lat,lon)

class Atmosphere:

    def __init__(self, year = 2000, doy = 1, model = 'advanced'):
        self.model = model
        self.year = year
        self.doy = doy
        self.air = Atmosphere.Air()
        self.rho_hist = []

    def at(self, timesec, altitude, latitude, longitude):
        if(self.model=='advanced'):
            output = atm.nrlmsise_output()
            inputatm = atm.nrlmsise_input(year = self.year, doy = self.doy, sec = timesec, alt=altitude*0.001, g_lat=latitude, g_long=longitude, lst=0.0, f107A=150., f107=150.0, ap=4.0, ap_a=None)
            flags = atm.nrlmsise_flags()
            flags.switches[:] = [1]*(len(flags.switches)) # output in kg/m³
            gtd7(inputatm,flags,output)
            rho = output.d[5]
            T = output.t[1]
        else:
            rho = density(altitude*0.001)
            T = 214.
        self.air.set_state(rho,T)
        self.rho_hist.append(rho)
        print('density at %f :', altitude, rho)
        return self.air

    class Air:

        def __init__(self):
            self.N_A = 6.023 * 10**23 # Avogadro
            self.MM_air = 0.0289644 # kg/mol
            self.sigma_c = 1e-19 # m² mean collisionnal cross section
            self.b_air = 1.458e-6 # kg/ms(K)^0.5
            self.S_air = 110.4 # K

        def set_state(self,rho,T):
            self.rho = rho # mass density
            self.nv = rho*self.N_A/self.MM_air # particle density
            self.mean_free_path = 1/(m.sqrt(2)*self.nv*self.sigma_c) # mean free path between collisions
            self.T = T # kinetic air temperature
            self.visc = self.b_air*self.T**(1.5)/(self.T+self.S_air) # Sutherland law for gases
            



class Projectile:
    
    def __init__(self, diameter, density):
        self.d = diameter
        self.S = m.pi*diameter**2
        self.m = m.pi/6.*diameter**3*density

    def propel(self, r_ini, v_ini):
        self.r_0 = r_ini
        self.r = r_ini.copy()
        self.v_0 = v_ini
        self.v = v_ini.copy()
        self.traj = np.array([r_ini])
        self.vel = np.array([v_ini])
        self.Kn = []

    def clocktick(self,dt,dv,lpm):
        self.v += dv
        self.r += self.v*dt
        self.Kn.append(lpm/self.d)
        self.traj = np.vstack((self.traj,self.r))
        self.vel = np.vstack((self.vel,self.v))

class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        print('called by :',cls)
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
 
class Mechanics(metaclass=Singleton):

    def __init__(self):
        self.animate = []
        print('new Mech :', self)
        self.clock = 0.
        self.timegens = []
        self.timeline = [self.clock]

    @classmethod
    def initialise(cls):
        cls._instances.clear()
        return Mechanics()

    @classmethod
    def getMech(cls):
        return Mechanics._instances[cls]

    def set(self,earth,schedule):
        self.earth = earth
        self.schedule = schedule

    def add_animate(self, obj):
        self.animate.append(obj)

    def fetch_events(self):
        evts = self.schedule.at(self.clock)
        for evt in evts:
            evt.trigger()
            print('Fetched : '+evt.name)
            tgen = evt.timegen
            print('Generator : ',tgen)
            self.timegens.append(tgen)

    def start(self):
        # 1 - list of events on schedule
        # 2 - calculation of time pace
        # 3 - update of animate objects
        self.fetch_events()
        # Loop until job finished
        while len(self.timegens)!=0:
            dts = []
            for i,tgen in enumerate(self.timegens.copy()):
                try:
                    dts.append(next(tgen))
                except StopIteration:
                    self.timegens.pop(i)
                    print('depleted')
                finally:
                    if(len(self.timegens)==0):
                        print('Mechanics stopped')
                        return
            dt = min(dts)
            self.move(dt)
            self.clock += dt
            self.timeline.append(self.clock)
            self.fetch_events()
                
    def move(self, dt):
        if(dt!=0):
            for obj in self.animate:
                if(isinstance(obj,Satellite)):
                    obj.clocktick(dt)
                elif(isinstance(obj,Projectile)):
                    dv, lpm = self.drag(obj,dt)
                    obj.clocktick(dt,dv,lpm)
                else:
                    obj.clocktick(dt)
    
    def N_coll(self,n_mean):
        if(n_mean >= 50):
            n_coll = int(np.random.normal(n_mean,0.5*n_mean))
        else:
            n_coll = 0
            for i in range(int(n_mean)*4):
                n_coll += int(rand.random()>=0.5)
        return n_coll

    def drag(self,obj,dt):
        volume = obj.S*np.linalg.norm(obj.v)*dt # Tube volume covered by projectile during dt
        n_volmean = self.earth.atmosphere(obj.r).nv
        # n_volmean = density((np.linalg.norm(obj.r)-self.earth.R)*0.001)*N_A/MM_air
        lpm = 1/(m.sqrt(2)*n_volmean*sigma_c) # mean free path formulae
        n_mean = n_volmean*volume
        n_coll = self.N_coll(n_mean)
        dv_coll = - n_coll*MM_air/(N_A*obj.m) * obj.v
        
        a_p = -obj.r/np.linalg.norm(obj.r)*mu_G/(np.linalg.norm(obj.r))**2

        dv_p = a_p*dt + dv_coll
        return dv_p, lpm


class Schedule:
    
    def __init__(self):
        self.events = {}
        self.times = []
    
    def plan(self, eventname, func, *args, time = 0, duration = iter([0])):
        evt = Event(eventname,duration,func,*args)
        if(time in self.events):
            self.events[time].append(evt)
        else:
            self.events[time] = [evt]
            self.times.append(time)
        self.times.sort()

    def at(self, time):
        if(len(self.times)!=0 and time>=self.times[0]):
            return self.events[self.times.pop(0)]
        else:
            return []

class Event:

    def __init__(self, name, timegen, func, *args):
        self.name = name
        self.func = func
        self.args = args
        self.timegen = timegen

    def trigger(self):
        print(self.name)
        if(self.func is not None):
            self.func(*self.args)   

class Reference:
    _refdate = dt.datetime(2000,1,1,12)

    @classmethod
    def get_refdate(cls):
        return Reference._refdate