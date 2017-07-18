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
from mathutils import Quaternion, Vector
import methutil as methu

global mu_G, R_G
mu_G = 398600.4418 * 10**9 # m³/s²
R_G = 6378.1 * 10**3 # m 


MM_air = 0.0289644 # kg/mol
N_A = 6.023 * 10**23
sigma_c = 1e-19 # m²

class Singleton(type):
    ''' Implementation of the Singleton design pattern on Python'''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        print('called by :',cls)
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        print(Singleton._instances)
        return cls._instances[cls]

class Earth(metaclass=Singleton):
    ''' Earth class encompasses the natural elements proper to earth :
    * Time frame(s) : Sidereal and Solar
    * Geometry,
    * Gravity field,
    * Atmosphere
    '''

    def __init__(self):
        self.R = R_G
        self.mu = mu_G
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
        self.targets = []
        self.lostraj = {}
        

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

class Projectile(Movody):
    ''' Passive object representing a solid sphere and recording its kinetic state trajectory '''

    smooth = 'smooth'
    coarse = 'coarse'

    def __init__(self, diameter, density, finition = smooth):
        self.d = diameter
        self.S = 0.25*m.pi*diameter**2
        self.m = m.pi/6.*diameter**3*density
        self.surftype = finition

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


class Orbit:
    ''' Classical keplerian orbit comprising :
    * semi-major axis,
    * eccentricity,
    * inclination,
    * ascending node argument,
    * perigea argument'''

    def __init__(self,params):
        self.a, self.e, self.i, self.Om, self.w = params
        self.n = m.sqrt(mu_G/self.a**3)
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
        p = np.sum(r_z**2)/mu_G
        a = 1/(2/r_0 - v_0**2/mu_G) # Vis-viva orbital energy equation
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
        v_x = -np.sqrt(mu_G/(self.a*(1-self.e**2)))*np.sin(theta)
        v_y = np.sqrt(mu_G/(self.a*(1-self.e**2)))*(self.e+np.cos(theta))
        v_z = np.array([0.]*v_x.size)
        return self.rel2abs.dot(np.array([v_x,v_y,v_z])).transpose()

    def getVel2(self,theta):
        ''' Get absolute velocity according to energy relation and tangent vector'''
        r = self.a*(1-self.e**2)/(1+self.e*m.cos(theta))
        vnorm = m.sqrt(mu_G*(2/r-1/self.a))
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
        self.air = Atmosphere.Air()
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
        self.air.set_state(rho,T)
        self.rho_hist.append(rho)
        return self.air

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

        def __init__(self):
            self.N_A = 6.023 * 10**23 # Avogadro
            self.MM_air = 0.0289644 # kg/mol
            self.R_m = 8.315/self.MM_air # Jmol/kgK
            self.gamma = 1.4 # gas caloric constant
            self.sigma_c = 1e-19 # m² mean collisionnal cross section
            self.b_air = 1.458e-6 # kg/m.s(K)^0.5
            self.S_air = 110.4 # K

        def set_state(self,rho,T):
            self.rho = rho # mass density
            self.nv = rho*self.N_A/self.MM_air # particle density
            self.mean_free_path = 1/(m.sqrt(2)*self.nv*self.sigma_c) # mean free path between collisions
            self.T = T # kinetic air temperature
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


    def grav(self,dt):
        a_p = self.earth.gravity(self.movody.r)
        # a_p2 = -self.earth.mu * self.movody.r/np.linalg.norm(self.movody.r)**3
        # print("Precise : ",a_p1, "\n", "Coarse : ", a_p2)
        return a_p*dt

    def N_coll(self,n_mean):
        n_coll = n_mean
        if(n_mean <= 10000):
            n_coll = int(np.round((np.random.normal(n_mean,1/m.sqrt(n_mean)))))
        return n_coll

    def drag(self,dt):
        volume = self.movody.S*np.linalg.norm(self.movody.v)*dt # Tube volume covered by projectile during dt
        
        air = self.earth.atmosphere(self.movody.r)
        n_volmean = air.nv
        # n_volmean = density((np.linalg.norm(obj.r)-self.earth.R)*0.001)*N_A/MM_air
        lpm = 1/(m.sqrt(2)*n_volmean*air.sigma_c) # mean free path formulae
        n_mean = n_volmean*volume
        n_coll = self.N_coll(n_mean)
        dv_coll = - n_coll*air.MM_air/(air.N_A*self.movody.m) * self.movody.v
        
        # a_p = -self.movody.r/np.linalg.norm(self.movody.r)*mu_G/(np.linalg.norm(self.movody.r))**2

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

class HypersonicBox(DynamicBox):
    ''' Helper class that supervises the dynamics of the projectile during its descent in the atmosphere'''

    _surftype = {Projectile.smooth: 0, Projectile.coarse : 1}

    def __init__(self, earth, time = 0):
        super().__init__(earth,time)
        self.sigma_d = 0
        self.machline = []
        self.rholine = []
        self.presline = []
        self.Knline = []
        self.Reline = []
        self.dragline = []

    def clocktick(self,dt):
        super().clocktick(dt)
        #dv_d, air_d = super().drag(dt)
        dv, air, C_D = self.drag(dt)
        dv_g = self.grav(dt)
        self.movody.clocktick(dt,dv+dv_g)
        self.rholine.append(air.rho)
        self.Knline.append(air.mean_free_path/self.movody.d)
        v_p = np.linalg.norm(self.movody.v)
        self.machline.append(v_p/air.c_a)
        # self.presline.append(air.p)
        self.Reline.append(air.rho*v_p*self.movody.d/air.visc)
        self.dragline.append(C_D)

    def set_movody(self,movody):
        super().set_movody(movody)
        self.sigma_d = HypersonicBox._surftype[movody.surftype]

    def dragsim(self,dt):
        volume = self.movody.S*np.linalg.norm(self.movody.v)*dt # Tube volume covered by projectile during dt

        air = self.earth.atmosphere(self.movody.r)
        n_volmean = air.nv
        
        # n_volmean = density((np.linalg.norm(obj.r)-self.earth.R)*0.001)*N_A/MM_air
        lpm = 1/(m.sqrt(2)*n_volmean*sigma_c) # mean free path formulae
        n_mean = n_volmean*volume
        n_coll = self.N_coll(n_mean)
        
        dv_coll = - n_coll*air.MM_air/(air.N_A*self.movody.m) * self.movody.v
    
        a_p = -self.movody.r/np.linalg.norm(self.movody.r)*mu_G/(np.linalg.norm(self.movody.r))**2

        dv_p = a_p*dt + dv_coll
        return dv_p, air

    def drag(self,dt):
        air = self.earth.atmosphere(self.movody.r)
        v_p = np.linalg.norm(self.movody.v)
        S = np.sqrt(air.gamma/2)*v_p/air.c_a
        C_D = m.e**(-S**2)/(np.sqrt(2)*S**3)*(1+2*S**2) + (2*S**4+2*S**2-1)/(2*S**4)*m.erf(S) + 2*self.sigma_d*np.sqrt(np.pi)/(3*S)
        dv_p = - 0.5*C_D*self.movody.S*air.rho*v_p*self.movody.v/self.movody.m*dt

        return dv_p, air, C_D

    # def N_coll(self,n_mean):
    #     if(n_mean >= 10):
    #         n_coll = int(np.random.normal(n_mean,1/m.sqrt(n_mean)))
    #     else:
    #         n_coll = 0
    #         for i in range(int(n_mean)*4):
    #             n_coll += int(rand.random()>=0.5)
    #     return n_coll


        
class TimeGenerators:

    def projectileTimeGenerator(projectile, N_itmax = 10**4):
        dt = 10**(-3)*R_G/np.linalg.norm(projectile.v_0)
        Kn = 10
        MM_air = 0.0289644
        n_it = 0
        alt = np.linalg.norm(projectile.r)-R_G

        while (Kn >= 1 and n_it < N_itmax and alt>=0):

            alt = np.linalg.norm(projectile.r)-R_G
            n_volmean = density(alt*0.001)*N_A/(MM_air)
            lpm = 1/(m.sqrt(2)*n_volmean*sigma_c)
            Kn = lpm/projectile.d
            n_it += 1

            yield dt

    def orbitalTimeGenerator(movody, duration, N_itmax = 10**4):
        dt = 10**(-3)*R_G/np.lin*alg.norm(projectile.v_0)
        t = 0

        while(t<duration and n_it < N_itmax):
            t += dt
            if(t>duration):
                dt -= t-duration

            yield dt
 
class Mechanics(metaclass=Singleton):
    ''' Master class for executing events, ticking the clock and notifying animated objects'''
    boxing = {Projectile: HypersonicBox}

    def __init__(self):
        self.animate = []
        print('new Mech :', self)
        self.clock = 0.
        self.timegens = []
        self.timeline = [self.clock]
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

    def set(self,earth,schedule):
        self.earth = earth
        self.schedule = schedule

    def add_animate(self, obj):
        if(obj.__class__ in Mechanics.boxing):
            box = Mechanics.boxing[obj.__class__](self.earth,self.clock)
            box.set_movody(obj)
            self.boxes[obj] = box   
            obj = box
            # self.hypersonicbox = HypersonicBox(self.earth,self.clock)
            # self.hypersonicbox.set_airobj(obj)
            # obj = self.hypersonicbox            
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
                    self.timegens.remove(tgen)
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
                    obj.clocktick(dt)


class Schedule:
    ''' Basic implementation of a schedule containing events disposed on a timetable. 
    Schedule can be consulted at any time for fetching events on schedule for one time only read.'''
    
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
    ''' Object containing a name, a task function and a slicing time generator'''

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
    ''' Convenient class for regrouping time and frame references of the code '''
    _refdate = dt.datetime(2000,1,1,12)

    @classmethod
    def get_refdate(cls):
        return Reference._refdate

    @classmethod
    def get_refloc(cls):
        return 35.6895*m.pi/180, 139.69171*m.pi/180, 'Tokyo' # Tokyo coordinates : Latitude, Longitude