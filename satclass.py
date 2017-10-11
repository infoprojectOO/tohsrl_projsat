# -*- coding: utf-8 -*-
"""
Created on Thu May 25 11:42:35 2017

@author: ASiapan
"""

import math as m
import param as cst
import numpy as np
import datetime as dt
from mathutils import Quaternion, Vector
import methutil as methu
from threading import Thread
import threading, traceback
import dsmc
from physenv import Earth, Atmosphere

class Singleton(type):
    ''' Implementation of the Singleton design pattern on Python'''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        # print('called by :',cls)
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        # print(Singleton._instances)
        return cls._instances[cls]

class Movody:
    '''Generic class for moving bodies = (movodies) in space'''

    def __new__(cls,*args,**kwargs):
        self = super().__new__(cls)
        self.specificities = {}
        return self

    def __init__(self, mass, surface, position0, velocity0, t0 = 0):
        ''' Create a new movody object with initial conditions and properties

        Arguments 
        ------
        mass : ( kg )
        surface : (m²)
        position0 : array[3] (m)
        velocity0 : array[3] (m/s)
        '''
        self.m = mass
        self.m0 = mass
        self.S = surface
        self.r = position0.copy()
        self.v = velocity0.copy()
        self.t = t0
        self.traj = np.array([position0]) 
        self.vel = np.array([velocity0])
        self.time = [self.t]

    def clocktick(self,dt,acc):
        ''' Generic function for updating any animated object through time
        
        Arguments 
        ------
        dt : seconds elapsed ( s )
        acc : the acceleration exerted on the body (m/s²)
        '''
        self.t += dt
        dv = acc*dt
        self.r += self.v*dt + dv*dt*0.5 #Euler forward integration
        self.v += dv
        self.traj = np.vstack((self.traj,self.r))
        self.vel = np.vstack((self.vel,self.v))
        self.time.append(self.t)

    def add_spec(self, spec, val):
        ''' Specify a special coded characteristic of this body (numeric, physical, etc.)'''
        self.specificities[spec] = val

class Thermovody(Movody):
    '''Generic class for thermodynamic bodies = (thermovodies) in space'''

    def __init__(self, mass, surface, temperature, *args):
        super().__init__(mass, surface, *args)
        self.T = temperature
        self.temp = [self.T]

    def clocktick(self,dt,acc,q):
        ''' Generic function for updating any thermally active object through time
        
        Arguments 
        ------
        dt : seconds elapsed ( s )
        acc : the acceleration exerted on the body (m/s²)
        q : the heat exerted on the body (J/s)
        '''
        super().clocktick(dt,acc)
        q_eff = q - cst.sigma_r*self.T**4*self.S*4
        self.T += q_eff*dt/(self.C_cal*self.m)
        try:
            assert self.T >= 0
        except(AssertionError):
            self.T = (q/(cst.sigma_r*self.S*4))**(0.25)
            print("Thermal saturation error !")
        self.temp.append(self.T)        

class Sphermovody(Thermovody):
    '''Generic class for spherical moving bodies = (sphermovodies) in space'''

    def __init__(self, diameter, material, *args):
        self.d = diameter
        self.material = material
        self.m_melt = 0
        self.rebulk()
        super().__init__(self.m, self.S, *args)
        self.diam = [self.d]

    def clocktick(self,dt,acc,q,Ni):
        ''' Generic function for updating any spherical object through time
        
        Arguments 
        ------
        dt : seconds elapsed ( s )
        acc : the acceleration exerted on the body (m/s²)
        q : the heat exerted on the body (J/s)
        Ni : the flux of impacting particle (#/s)
        '''
        super().clocktick(dt,acc,q)
        if self.ablation != 0: # Check if object is ablatable or not
            self.ablation = self.compute_ablation_rate(self.T,self.material.melting_T)
            n_abl = self.ablation*Ni*dt
            m_abl = n_abl*self.material.mol_mass/cst.N_A
            if(self.T >= self.material.melting_T):
                DT = self.T-self.material.melting_T
                self.m_melt += DT*(self.C_cal*self.m)/self.material.melt_specheat
                if(self.m_melt<=self.m0):
                    self.T = self.material.melting_T
                # else:
                    # m_abl += m_melt
            V_abl = m_abl/self.material.mass_density
            dV = self.V-V_abl
            if dV >= 0 : self.d = 2*(dV*3/(4*m.pi))**(1/3)
            else: self.d = 0
            self.rebulk()
        self.diam.append(self.d)

    def rebulk(self):
        ''' Update surface, volume and mass for any change in size '''
        self.S = 0.25*m.pi*self.d**2
        self.V = 0.5/3*m.pi*self.d**3
        self.m = m.pi/6.*self.d**3*self.material.mass_density

    def compute_ablation_rate(self,T,T_melt):
        ''' Set the ratio of particles ablated according to surface temperature and melting temperature

        For every impacting particle on the surface, regardless of its position and velocity, 
        a certain number of fraction of particles is emitted from the surface into the flow

        Arguments
        ---------
        T : body surface temperature (K)
        T_melt : body material melting point (K)
        '''
        return 10*m.e**((T-T_melt)/T_melt)

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
        self.camera = Camera(width*0.8) # The camera is automatically placed on the z- face of the satellite
        self.att0 = Quaternion([0.,1.,0.],m.pi*0.5) # Initial attitude of the satellite in the satellite reference frame
        self.setattitude()
        self.attitudes = [self.attitude.copy()]
        self.targets = []
        self.lostraj = {}
        self.clock = 0.        

    def clocktick(self, dt):
        ''' Generic function for animated objects moved accross dt (s) lapse'''
        self.clock += dt
        self.M += self.orbit.n*dt
        self.M = self.M%(2*m.pi)
        self.nu = self.orbit.anomaly(self.M,'true')
        self.setattitude()
        if(self.tracking):
            self.track()

    def eject(self, projectile, dv):
        ''' Eject the projectile according to the ejection velocity in the satellite-orbital frame

        Arguments 
        ---------
        projectile : Projectile object - the one being ejected
        dv : array[3] - (norm, theta, phi) velocity from the satellite in the orbital Fresnel polar coordinates'''
        (norm,theta,phi) = dv
        dvx = norm*(m.cos(theta)) # x points towards the satellite's direction
        dvy = norm*(m.sin(theta)*m.cos(phi)) # y points inwards in the orbital plane
        dvz = norm*(m.sin(theta)*m.sin(phi))  # z completes the triorthon
        r = self.orbit.getPos(self.nu)
        v = self.orbit.getVel(self.nu)+self.orbit.rel2abs.dot(self.orbit.getVel2Rel(self.nu).dot(np.array([dvx,dvy,dvz]).T)).T
        projectile.propel(r,v,self.clock)
        self.addTarget(projectile)
        # Mechanics.getMech().add_animate(projectile)
        return projectile

    def addTarget(self, projectile):
        ''' Adds a projectile on the list of tracked objects '''
        self.targets.append(projectile)
        self.tracking = True
        r = self.orbit.getPos(self.nu)
        self.linesight = projectile.r-r
        if(np.linalg.norm(self.linesight)!=0.):
            self.linesight = self.linesight/np.linalg.norm(self.linesight)
        self.lostraj[projectile] = np.array([self.linesight.copy()]) # Pointing direction in the absolute frame
        self.traj = np.array([r.copy()])
        #self.setattitude(-self.linesight)

    def track(self):
        ''' Update all the pointing directions towards the tracked objects '''
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
        '''Update satellite attitude in absolute coordinates 

        All rotations are rotated according to Quaternion assembly
        '''
        self.attitude = self.att0.copy()
        Om_rot = Quaternion([0.,0.,1.],self.orbit.Om) # Precession
        i_rot = Quaternion([1.,0.,0.],self.orbit.i) # Nutation
        i_rot.rotate(Om_rot) # Update the rotation axis of nutation according to precession changes
        wnu_rot = Quaternion([0.,0.,1.],self.orbit.w+self.nu) # Rotation
        wnu_rot.rotate(i_rot) # Update the rotation axis of rotation according the nutation & precession
        self.attitude.rotate(wnu_rot) # Update the satellite attitude according to the orbital position related rotation
        # pointing = np.array(orbit.getPos(self.nu)/np.linalg.norm(orbit.getPos(self.nu)))
        # rot_vec = np.cross(ref,pointing)
        # rot_norm = np.linalg.norm(rot_vec)
        # self.attitude = Quaternion(rot_vec/rot_norm,m.asin(rot_norm))

class Camera:

    def __init__(self,width):
        self.width = width

class Projectile(Sphermovody):
    ''' Passive object representing a solid sphere and recording its kinetic and internal state trajectory '''

    smooth = 'smooth'
    coarse = 'coarse'

    def __init__(self, diameter, material, T_surf = 200, finition = coarse, ablation = 0, id = None):
        self.d = diameter
        self.material = material
        self.surftype = finition
        self.ablation = ablation
        self.C_cal = material.caloric_capacity
        self.k_c = material.caloric_conductivity
        self.T = T_surf
        self.specificities = {DynamicBox.DSMC:True}
        self._id = id

    def propel(self, r_ini, v_ini, t_ini = 0):
        ''' Adopt position and velocity at an initial time in space '''
        super().__init__(self.d,self.material,self.T,r_ini,v_ini,t_ini)
        self.r_0 = r_ini.copy()
        self.v_0 = v_ini.copy()

    def clocktick(self,dt,acc,q,Ni):
        super().clocktick(dt,acc,q,Ni)

    @property
    def id(self): # Used for distinguishing projectiles
        return str(id(self))

    @property
    def id_given(self):
        ''' Used for naming projectiles at their creation '''
        return (self._id != None)

    def __str__(self):
        if(self._id): return self._id
        if(self.ablation):
            return 'Projectile with ablation'
        else:
            return 'Projectile with durability'

    def copy(self):
        return Projectile(self.d,self.material,self.T,self.surftype,self.ablation)

class Ghostectile(Sphermovody):
    ''' Copy version of a Projectile, foresees the trajectory in the transitional regime'''

    def __init__(self,projectile):
        self.real = projectile
        p = projectile
        super().__init__(p.d, p.material, p.T, p.r, p.v, p.t)
        self.surftype = p.surftype
        self.C_cal = p.C_cal
        self.k_c = p.k_c
        self.d = p.d
        self.r_0 = p.r.copy()
        self.v_0 = p.v.copy()
        self.ablation = p.ablation

    def copy(self):
        return Ghostectile(self)

    @property
    def id(self):
        return str(id(self))

    def __str__(self):
        return 'Ghost with approximate data'


class DynamicBox:
    ''' Helper class based on the Decorator Design Pattern for implementing dynamics on objects'''

    OVERRIDE_DRAG = 'drag'
    OVERRIDE_AERODOM = 'aerodom'
    OVERRIDE_KN_LIM = 'knudsen threshold'
    DSMC = 'compute with dsmc'

    def __init__(self, earth, time = 0):
        self.earth = earth
        self.atmosphere = earth.atm
        self.timeline = [time]
        self.clock = time

    def set_movody(self,movody):
        ''' Set the movody that this Box wraps into an aerodynamic regime'''
        self.movody = movody
        self.altline = np.array([(np.linalg.norm(self.movody.r)-Earth.R)*0.001])
        self.speedline = np.array([np.linalg.norm(self.movody.v)])

    def clocktick(self,dt):
        ''' Update on the motion and state of the movody '''
        self.clock += dt
        self.timeline.append(self.clock)
        assert self.movody is not None
        self.altline = np.hstack((self.altline,(np.linalg.norm(self.movody.r)-Earth.R)*0.001))
        self.speedline = np.hstack((self.speedline,np.array([np.linalg.norm(self.movody.v)])))

    def index_at(self,alt):
        ''' Return the index of the array recordings corresponding to the altitude given

        Arguments
        ---------
        alt : altitude (km)

        Returns
        -------
        int - index
        '''
        if(type(alt) is list):
            inds = []
            for a in alt:
                inds.append(methu.find_closest(self.altline,a))
            return inds
        else:
            return methu.find_closest(self.altline,alt)

    def grav(self):
        ''' Compute the gravity acceleration on the movody'''
        a_p = self.earth.gravity(self.movody.r)
        # a_p2 = -self.earth.mu * self.movody.r/np.linalg.norm(self.movody.r)**3
        # print("Precise : ",a_p1, "\n", "Coarse : ", a_p2)
        return a_p

    def N_coll(self,n_mean):
        ''' Generate a randomly chosen number of particle collisions '''
        n_coll = n_mean
        if(n_mean <= 1000):
            n_coll = int(np.round((np.random.normal(n_mean,1/m.sqrt(n_mean)))))
        return n_coll

    def drag(self,dt):
        ''' Compute the drag exerted on the movody'''
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

    def set_movody(self,satellite):
        self.sat = satellite


class HypersonicBox(DynamicBox):
    ''' Helper class that supervises the dynamics of the projectile during its descent in the atmosphere'''

    def __init__(self, earth, time = 0):
        ''' Initialise the Box with all the necessary arrays extended to record timely data during the projectile descent'''
        super().__init__(earth,time)
        self.Kn_trans = 10
        self.machline = []
        self.rholine = []
        self.presline = []
        self.trline = []
        self.Knline = []
        self.Reline = []
        self.dragline = []
        self.heatline = []
        self.collline = []
        self.regimes = {}
        self.ghostbox = None
        self.ghost = None

    def clocktick(self,dt):
        if(self.movody.d == 0): 
            Mechanics().remove_animate(self)
            return
        super().clocktick(dt)
        #dv_d, air_d = super().drag(dt)
        # Compute the drag and gravity accelerations
        air = self.earth.atmosphere(self.movody.r)
        self.regime = AeroRegime(self.movody,air)
        a_d, C_D, q, Ni = self.drag(self.regime)
        a_g = self.grav()
        
        # Record the new regime values
        self.rholine.append(air.rho)
        self.Knline.append(self.regime.Kn_mix)
        self.machline.append(self.regime.Ma)
        self.presline.append(air.p)
        self.trline.append(self.regime.tr)
        self.Reline.append(air.rho*self.regime.U*self.movody.d/air.visc)
        self.dragline.append(C_D)
        self.heatline.append(q)
        self.collline.append(Ni)

        # Update the movody
        self.movody.clocktick(dt,a_d+a_g,q,Ni)
        

    def set_movody(self,movody):
        ''' Set the proper aerodynamic domains scheme according to the movody's specificities'''
        super().set_movody(movody)
        self.Kn_trans = movody.specificities.get(self.OVERRIDE_KN_LIM,self.Kn_trans)
        # self.ablate = movody.specificities.get(self.ABLATION,0)
        if(self.OVERRIDE_DRAG in movody.specificities): # Drag coefficient is constant and overrided by a given value
            self.regimes[AeroDom.Custom] = AeroDom(movody.specificities[self.OVERRIDE_DRAG])
            self.movody.aeroreg = AeroDom.Custom
        elif(self.OVERRIDE_AERODOM in movody.specificities): # The Aerodynamic Domain is overrided by only one for the whole descent
            self.regimes[AeroDom.Custom] = AeroDom(movody.specificities[self.OVERRIDE_AERODOM])
            self.movody.aeroreg = self.regimes[AeroDom.Custom].name
        elif(self.DSMC in movody.specificities): # DSMC is used for the movody simulation
            movody.aeroreg = AeroDom.Free_Molecular_Flow
            if(movody.id_given): id = str(movody)
            else: id = None
            self.dsmc = dsmc.DSMC(id = id)
            self.regimes[AeroDom.Free_Molecular_Flow] = FMF()
            self.regimes[AeroDom.Transitional_Flow] = TF(self.dsmc)
            self.regimes[AeroDom.Near_Continuum_Flow] = self.regimes[AeroDom.Transitional_Flow]
        else: # Only theoretical models are used for the descent
            movody.aeroreg = AeroDom.Free_Molecular_Flow
            self.regimes[AeroDom.Free_Molecular_Flow] = FMF()
            self.regimes[AeroDom.Transitional_Flow] = TF()
            self.regimes[AeroDom.Near_Continuum_Flow] = NCF()

    def drag(self,regime):
        ''' Compute the drag of the aerodynamic regime'''
        dom = self.find_dom(regime) # find the domain (FMF/ TF/ NCF ?)
        C_D = dom.get_dragC(regime) # retrieve the drag, 
        q = dom.get_heat(regime)  # heat
        Ni = dom.get_onflow(regime) # and impact flux
        acc_drag = - 0.5*C_D*self.movody.S*regime.air.rho*regime.U*self.movody.v/self.movody.m # deduce the deceleration caused by drag
        return acc_drag, C_D, q, Ni

    def find_dom(self,regime):
        ''' Find the appropriate domain corresponding to an aerodynamic regime instance'''
        if(self.OVERRIDE_DRAG in self.movody.specificities) : return self.regimes[AeroDom.Custom]
        if(self.OVERRIDE_AERODOM in self.movody.specificities) : return self.regimes[AeroDom.Custom]
        Kn_mix = regime.Kn_mix

        # Compute hypersonic flow compression ratios (estimates)
        Ma2 = regime.Ma**2 
        gamma = regime.air.gamma
        compression_rho = (gamma+1)*Ma2/((gamma-1)*Ma2 + 2)

        Kn_s = regime.Kn/compression_rho # Knudsen after the Bow shock
        if self.DSMC in self.movody.specificities and Kn_mix <= self.Kn_trans+10 and (self.ghostbox is None): # Ghost used only for DSMC cases, when the mixed Knudsen value reaches below the threshold
            self.ghost = Ghostectile(self.movody)
            self.namecode = 'Ghost Projectile Cloned'
            Mechanics().schedule.urge(self.namecode, self.ghostready, delay = 0., duration = TimeGenerators().ghostectileTimeGenerator(self.ghost),type = Mechanics.EVENT_PRIOR, key = self.ghost.id)
        if(Kn_mix >= self.Kn_trans): # Free molecular flow case
            self.movody.aeroreg = AeroDom.Free_Molecular_Flow
            return self.regimes[AeroDom.Free_Molecular_Flow]            
        elif(Kn_s>=0.05): # Transitional flow case
            self.movody.aeroreg = AeroDom.Transitional_Flow
            return self.regimes[AeroDom.Transitional_Flow]
        else: # Near continuum flow case
            self.movody.aeroreg = AeroDom.Near_Continuum_Flow
            return self.regimes[AeroDom.Near_Continuum_Flow]

    def ghostready(self):
        ''' Launch ghost '''
        self.ghostbox = Mechanics().add_animate(self.ghost,self.ghost.id)
        self.ghostbox.set(self.dsmc)

class GhostsonicBox(HypersonicBox):
    ''' Manage the descent of the ghostectile and launch DSMC simulations'''

    def __init__(self, earth, time = 0):
        super().__init__(earth,time)
        self.pathmark = []
        self.markindex = []

    def clocktick(self,dt):
        super().clocktick(dt)
        if self.differs(): # New DSMC should be launched
            self.record(self.regime) # record the input parameters used for DSMC
            self.request_sim()

    def differs(self):
        ''' Decide whether a new DSMC simulation should be launched based on the previous regime instance
        
        Checks if the current regime differs from the previous to a relative degree of ...%

        Returns
        -------
        boolean <- Should a new DSMC be launched ?
        '''
        prev = self.pathmark[-1]
        now = self.regime
        # alt_now = (np.linalg.norm(self.movody.r)-self.earth.R)*0.001
        # alt_prev = self.pathmark[-1].alt
        return prev.rel_dist(now) >= 0.1 # Aeroregime change of at least 10 %

    def record(self,aeroreg):
        ''' Record the regime as an instance for DSMC input database'''
        self.markindex.append(max(len(self.rholine)-1,0))
        self.pathmark.append(aeroreg)

    def request_sim(self):
        self.dsmc.add_sim(self.pathmark[-1],dsmc.SimInput.TIMESAMP.default)
        # pass

    def set_movody(self,movody):
        super().set_movody(movody)
        self.record(AeroRegime(self.movody,self.earth.atmosphere(self.movody.r)))

    def set(self,dsmc):
        self.dsmc = dsmc
        self.request_sim()

class AeroRegime:
    ''' Represents the aerothermodynamic regime : interaction of air with hypersonic thermovody '''

    t_w = 0.1
    Ma_w = 0.3
    Re_w = 0.6
    Kn_w = 0.
    surftype = {Projectile.smooth: 0, Projectile.coarse : 1}

    def __init__(self,thermovody,air):
        self.air = air
        self.alt = (np.linalg.norm(thermovody.r)-Earth.R)*0.001
        self.U = np.linalg.norm(thermovody.v)
        self.L = thermovody.d # Characteristic object size
        self.m = thermovody.m
        self.C_cal = thermovody.C_cal
        self.S = thermovody.S # Characteristic object section area
        self.Tr = thermovody.T # Reflected flow temperature
        self.tr = self.Tr/self.T # temperature ratio (reflected/oncoming)
        self.sigma_d = self.surftype.get(thermovody.surftype,thermovody.surftype) # surface accommodation coefficient [0;1]
        self.Ma = self.U/self.air.c_a 
        self.T0 = air.T*(1+0.5*(air.gamma-1)*self.Ma**2) # Stagnation point temperature (total)
        self.Re = self.U*self.L*self.air.rho/self.air.visc
        self.v_diff = np.sqrt(8*air.R_m*thermovody.T/m.pi) # mean diffused particle velocity
        self.Kn = self.air.mean_free_path/self.L
        self.Kn_mix = self.Kn*self.v_diff*np.sqrt(2)/self.U # lambda_m = v_moy/(n*U*sigma)

    @property
    def T(self):
        return self.air.T
    @property
    def nv(self):
        return self.air.nv
    @property
    def p(self):
        return self.air.p

    @property
    def Re_0(self):
        return self.Re*np.sqrt(self.air.T/self.T0)

    @property
    def q(self):
        return 0.5*self.C_H*self.air.rho*self.S*self.U**3
        #self.xq*self.S*self.p*self.air.v_therm
    @property
    def Ni(self):
        return 0.5*self.C_Ni*self.nv*self.S*self.U
        #self.xNi*self.S*self.nv*self.air.v_therm

    def norm(self):
        ''' Compute the "norm" of this regime based on similarity parameters'''
        return np.sqrt(self.Re**2*self.Re_w + self.Ma**2*self.Ma_w + self.tr**2*self.t_w + self.Kn**2*self.Kn_w)

    def distance(self,reg):
        ''' Compute the inter-regime distance according to similarity parameters

        Arguments
        ---------
        reg : AeroRegime object - other regime to differentiate from
        '''
        return np.sqrt((self.Re-reg.Re)**2*self.Re_w + (self.Ma-reg.Ma)**2*self.Ma_w + (self.tr-reg.tr)**2*self.t_w + (self.Kn-reg.Kn)**2*self.Kn_w)

    def rel_dist(self,reg):
        ''' Relative interregime distance'''
        return self.distance(reg)/self.norm()

    def copy(self):
        thermovody = Thermovody(self.m,self.S,self.Tr,[self.alt*1000+Earth.R,0,0],[self.U,0,0])
        thermovody.d = self.L
        thermovody.C_cal = self.C_cal
        thermovody.surftype = self.sigma_d
        return AeroRegime(thermovody,self.air.copy())


class AeroRegimeFactory:
    ''' Factory design pattern class for easily building up regime objects from various input parameters'''

    def __init__(self,earth):
        self.atm = earth.atm
        self._Ma = None
        self._Re = None
        self._Kn = None
        self.air = None
        self.T = None
        self.set_default()

    def set_default(self, alt=100, speed=7000):
        self.projectile = Projectile(cst.d_p,cst.Material.COPPER,T_surf=300)
        self.projectile.propel([alt*1000+Earth.R,0,0],[speed,0,0])

    def at(self,altitude):
        ''' Set the altitude to get the air conditions 

        Arguments
        ---------
        altitude (km)
        '''
        self.air = self.atm.at(0.,altitude*1000,0,0)

    def with_(self,projectile):
        ''' Give a projectile instance into the regime '''
        self.projectile = projectile
        
    @property
    def Ma(self):
        return self._Ma
    @property
    def Re(self):
        return self._Re
    @property
    def Kn(self):
        return self._Kn

    @Ma.setter
    def Ma(self,Ma):
        self._Ma = Ma

    @Re.setter
    def Re(self,Re):
        self._Re = Re

    @Kn.setter
    def Kn(self,Kn):
        self._Kn = Kn

    def create(self):
        ''' Create new regime object with the input parameters given to the factory object '''
        assert self.projectile is not None # Cannot operate without a projectile instance
        L = self.projectile.d
        air = self.air
        if air is None: # Build up air's state from similarity parameters 
            T = self.T
            if(not T): T = self.projectile.T
            c_a = np.sqrt(cst.gamma*cst.R_m*T)
            v_therm = np.sqrt(2*cst.R_m*T)
            visc = cst.b_air*T**(1.5)/(T+cst.S_air)
            # assert not (self._Ma and self._Kn and self._Re)
            if self._Re and self._Ma and self._Kn:
                pass
            if self._Re and self._Ma:
                U = self._Ma*c_a
                rho = self._Re*visc/(U*L)
            elif self._Kn and self._Ma:
                U = self._Ma*c_a
                n = v_therm/(U*cst.sigma_c*self._Kn)
                rho = cst.MM_air*n/cst.N_A
            elif self.Kn and self.Re:
                raise(Exception('Not implemented for Kn and Re together !'))
            else:
                raise(Exception('Not enough input arguments'))
            air = Atmosphere.Air(rho,T)
            self.projectile.propel(np.array([100000,0,0]),np.array([U,0,0]))
        return AeroRegime(self.projectile,air)

# Strategy Desgin pattern for dividing the different flow regimes in the box.
class AeroDom:
    ''' Classes distinguishing the three rarefied gas hypersonic domains'''

    Free_Molecular_Flow = 'fmf'
    Transitional_Flow = 'tf'
    Near_Continuum_Flow = 'ncf'
    Custom = "?"

    def __init__(self, overC_D = 2):
        self.C_D = overC_D

    def get_dragC(self,aeroreg):
        '''Drag coefficient C_D '''
        return self.C_D

    def get_heat(self,aeroreg):
        ''' Heat flow integrated over the whole body '''
        return 0

    def get_onflow(self,aeroreg):
        ''' Oncoming particle flux over the whole body '''
        U = aeroreg.U
        N_i = aeroreg.S*aeroreg.air.nv*U
        return N_i

class FMF(AeroDom):
    ''' Free Molecular Flow theoretical implementation'''

    @property
    def name(self):
        return AeroDom.Free_Molecular_Flow

    def get_dragC(self,aeroreg):
        U = aeroreg.U
        S = np.sqrt(aeroreg.air.gamma/2)*aeroreg.Ma
        Sr = U/np.sqrt(2*aeroreg.air.R_m*aeroreg.Tr)
        sigma_d = aeroreg.sigma_d
        C_D = m.e**(-S**2)/(np.sqrt(m.pi)*S**3)*(1+2*S**2) + (4*S**4+4*S**2-1)/(2*S**4)*m.erf(S) + 2*sigma_d*np.sqrt(np.pi)/(3*Sr)
        return C_D

    def get_heat(self,aeroreg):
        U = aeroreg.U
        S = np.sqrt(aeroreg.air.gamma/2)*aeroreg.Ma
        p = aeroreg.air.p # rho R T
        sqre = np.sqrt(aeroreg.air.R_m*aeroreg.T/2) # n* sqrt(kT/2 m)
        diff = aeroreg.sigma_d
        Tr = aeroreg.Tr
        zeta = 4
        q = aeroreg.S*diff*p*sqre* ((2*S**2+1+(4+zeta)*(1-Tr/aeroreg.T))*(m.erf(S)/S*(S**2+0.5)+m.e**(-S**2)/(m.sqrt(m.pi)))-m.erf(S)/S)
        return q

    def get_onflow(self,aeroreg):
        U = aeroreg.U
        S = np.sqrt(aeroreg.air.gamma/2)*aeroreg.Ma
        v_therm = aeroreg.air.v_therm*np.sqrt(m.pi)/2
        N_i = aeroreg.S*aeroreg.air.nv*v_therm*((S**2+0.5)*m.erf(S)/S + m.exp(-S**2)/m.sqrt(m.pi))
        return N_i

class TF(AeroDom):
    ''' Transitional Flow semi-empirical and DSMC implementation'''

    @property
    def name(self):
        return AeroDom.Transitional_Flow

    def __init__(self,dsmc=None):
        self.dsmc = dsmc
        if not dsmc:
            self.fmf = FMF()
            self.ncf = NCF()
        self.dragCall = False
        self.heatCall = False

    def get_dragC(self,aeroreg):
        self.dragCall = not self.dragCall
        if(self.dragCall != self.heatCall):
            if(self.dsmc): self.extract(aeroreg)
            else: self.bridge(aeroreg)            
        return self.C_D

    def get_heat(self,aeroreg):
        self.heatCall = not self.heatCall
        if(self.dragCall != self.heatCall):
            if(self.dsmc): self.extract(aeroreg)
            else: self.bridge(aeroreg)
        return self.q

    def get_onflow(self,aeroreg):
        return self.N_i

    def bridge(self,aeroreg):
        Ma2 = aeroreg.Ma**2
        gamma = aeroreg.air.gamma
        compression_rho = (gamma+1)*Ma2/((gamma-1)*Ma2 + 2)
        rho_s = compression_rho*aeroreg.air.rho
        T_s = (2*gamma*Ma2-(gamma-1))*((gamma-1)*Ma2+2)/((gamma+1)**2*Ma2)*aeroreg.T
        if(T_s)>10000: T_s = 10000
        air_s = Atmosphere.Air(rho_s,T_s)
        factor = 0.1
        survive_rate = m.exp(-factor*compression_rho*aeroreg.L/air_s.mean_free_path)
        C_Dfmf = self.fmf.get_dragC(aeroreg)
        heatfmf = self.fmf.get_heat(aeroreg)
        Nifmf = self.fmf.get_onflow(aeroreg)
        C_Dnc = self.ncf.get_dragC(aeroreg)
        heatnc = self.ncf.get_heat(aeroreg)
        Ninc = self.ncf.get_onflow(aeroreg)

        self.C_D = C_Dnc + (C_Dfmf-C_Dnc)*survive_rate
        self.q = heatnc + (heatfmf-heatnc)*survive_rate
        self.N_i = Ninc + (Nifmf-Ninc)*survive_rate

    def extract(self,aeroreg):
        self.C_D, self.C_H, self.C_Ni = self.dsmc.exertion(aeroreg)
        self.q = self.C_H*0.5*aeroreg.S*aeroreg.air.rho*aeroreg.U**3
        self.N_i = self.C_Ni*0.5*aeroreg.S*aeroreg.nv*aeroreg.U


class NCF(AeroDom):
    ''' Near Continuum Flow theoretical approximate relations '''

    @property
    def name(self):
        return AeroDom.Near_Continuum_Flow

    def get_dragC(self,aeroreg):
        return 1

    def get_heat(self,aeroreg):
        Ma2 = aeroreg.Ma**2
        gamma = aeroreg.air.gamma
        T_s = (2*gamma*Ma2-(gamma-1))*((gamma-1)*Ma2+2)/((gamma+1)**2*Ma2)*aeroreg.T        
        if(T_s)>8000: T_s = 8000
        rho_s = aeroreg.air.rho*(gamma+1)*Ma2/((gamma-1)*Ma2 + 2)
        S = ((gamma-1)*Ma2 + 2)/(2*gamma*Ma2-(gamma-1))*np.sqrt(aeroreg.T*gamma/(T_s*2))
        p = rho_s*T_s*aeroreg.air.R_m # rho R T
        sqre = np.sqrt(aeroreg.air.R_m*T_s/2) # n* sqrt(kT/2 m)
        diff = aeroreg.sigma_d
        Tr = aeroreg.Tr
        zeta = 4
        q = aeroreg.S*diff*p*sqre*((2*S**2+1+(4+zeta)*(1-Tr/T_s))*(m.erf(S)/S*(S**2+0.5)+m.e**(-S**2)/(m.sqrt(m.pi)))-m.erf(S)/S)
        return q

    def get_onflow(self,aeroreg):
        Ma2 = aeroreg.Ma**2
        gamma = aeroreg.air.gamma
        T_s = (2*gamma*Ma2-(gamma-1))*((gamma-1)*Ma2+2)/((gamma+1)**2*Ma2)*aeroreg.T
        if(T_s)>8000: T_s = 8000
        rho_s = aeroreg.air.rho*(gamma+1)*Ma2/((gamma-1)*Ma2 + 2)
        S = ((gamma-1)*Ma2 + 2)/(2*gamma*Ma2-(gamma-1))*np.sqrt(aeroreg.T*gamma/(T_s*2))
        v_therm = aeroreg.air.v_therm*np.sqrt(T_s/aeroreg.T)*np.sqrt(m.pi)/2
        N_i = aeroreg.S*(rho_s/aeroreg.air.MM_air*aeroreg.air.N_A)*v_therm*((S**2+0.5)*m.erf(S)/S + m.exp(-S**2)/m.sqrt(m.pi))
        return N_i

class TimeGenerators(metaclass=Singleton):
    ''' Functions are generators which yield a lapse dt until flight conditions apply no longer'''

    def __init__(self):
        self.earth = Earth()
        self.atm = self.earth.atm

    def get_dt(self,projectile):
        ''' Choose a lapse according to aerodynamic domain '''
        if(projectile.aeroreg==AeroDom.Free_Molecular_Flow):
            dt = 10**(-3)*self.earth.R/np.linalg.norm(projectile.v)
        elif(projectile.aeroreg== AeroDom.Transitional_Flow):
            air = self.earth.atmosphere(projectile.r)
            v_p = np.linalg.norm(projectile.v)
            t_coll = 1/(air.nv*v_p*air.sigma_c)
            dt = t_coll*1000
        else:
            dt = 10**(-3)*self.earth.R/np.linalg.norm(projectile.v)
        return dt

    def ghostectileTimeGenerator(self,ghostectile, N_itmax = 10**4):
        alt = (np.linalg.norm(ghostectile.r)-self.earth.R)*0.001
        n_it = 0
        while (alt >= 50 and ghostectile.d>0 and n_it < N_itmax and alt>=0):
            alt = (np.linalg.norm(ghostectile.r)-self.earth.R)*0.001
            dt = 10**(-3)*self.earth.R/np.linalg.norm(ghostectile.v)
            n_it += 1
            yield dt

    def projectileTimeGenerator(self,projectile, N_itmax = 10**5):
        alt = (np.linalg.norm(projectile.r)-self.earth.R)*0.001
        n_it = 0
        while (alt >= 50 and projectile.d > 0 and n_it < N_itmax and alt>=0):
            alt = (np.linalg.norm(projectile.r)-self.earth.R)*0.001
            dt = self.get_dt(projectile)
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
    ''' Master class for executing events, ticking the clock and notifying animated objects

    The class allows for events to unfold on multiple threads, that have separate processing.
    All ordinary events will however be placed onto the mainline which decides when the mechanics stop or not'''
    boxing = {Projectile: HypersonicBox, Ghostectile: GhostsonicBox}
    EVENT_SEVERED = 'severed'
    EVENT_PRIOR = 'prior'

    def __init__(self):
        self.animate = []
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
        ''' Adds an animated object onto the mechanics update list

        Arguments
        ---------
        obj : Object with clocktick - receives clock ticks from the mechanics
        ---------- Keyword args
        key : hashable object that relates to a special thread
        '''
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

    def remove_animate(self,obj,key=0):
        ''' Adds an animated object onto the mechanics update list

        Arguments
        ---------
        obj : Object to be removed
        ---------- Keyword args
        key : hashable object that relates to a special thread from which the object is updated
        '''
        try:
            self.animate.remove(obj)
            self.timelines[self.linekeys[key]].animate.remove(obj)
        except:
            traceback.print_exc()

    def fetch_events(self):
        ''' Check the schedule if there are new events on the main agenda, prepare their timeline and trigger them'''
        evts = self.schedule.at(self.clock) # events due now
        for evt in evts:
            tgen = evt.timegen
            if(evt.type == self.EVENT_SEVERED or evt.type == self.EVENT_PRIOR): # Prepare a new timeline for the event
                self.linekeys[evt.key] = len(self.timelines)
                newtl = Mechanics.Timeline(self.clock)
                newtl.timegens.append(tgen)
                self.timelines.append(newtl) 
            evt.trigger()
            if(evt.type == self.EVENT_SEVERED): # Start the timeline independently from main            
                newtl.start()
            elif(evt.type == self.EVENT_PRIOR): # Ask main timeline to pause and join after the new timeline finishes
                newtl.start()
                self.timelines[self.linekeys[0]].affix(newtl)
            else: # Added to main timeline
                self.mainline.timegens.append(tgen)
                self.timegens.append(tgen)

    def notify(self, dt):
        '''Ask mechanics to look up for new events and update its clock time '''
        self.clock += dt
        self.fetch_events()
        return 'ok'

    def alert(self):
        ''' Warn the mechanics to leap the schedule onto the next task

        Returns
        -------
        boolean <- does the schedule have an upcoming event ?'''
        time = self.schedule.next()
        if time: # Time leap
            self.clock = time
            self.fetch_events()
        return (time != None)          

    def start(self):
        ''' Starts the main timeline and consults the schedule where to start '''
        self.fetch_events()
        self.mainline.run()

    def eject(self,sat,proj,vel):
        ''' Ask the satellite to eject a projectile

        Arguments
        ---------
        sat : Satellite object - satellite ALE
        proj : Projectile object - to be ejected next
        vel : array[3] - velocity in polar orbital coordinates (m/s)
        '''
        sat.eject(proj,vel)
        if(proj.id in self.linekeys):
            self.add_animate(proj,proj.id)
        else:
            self.add_animate(proj)


    class Timeline(Thread):
        ''' Represents a clock that records its time history and calls on to its animated objects'''

        def __init__(self,clockset = 0, call = None):
            super().__init__()
            self.clock = clockset
            self.timeline = [self.clock]
            self.animate = []
            self.timegens = []
            self.history = {}
            self.call = call

        def run(self):
            # 1 - list of events on schedule
            # 2 - calculation of time pace
            # 3 - update animate objects
            # Loop until job finished
            while len(self.timegens)!=0: # Still an ongoing time lapse
                loop = False
                dts = [] # Concurrent time lapses
                for i,tgen in enumerate(self.timegens.copy()):
                    try:
                        dts.append(next(tgen))
                    except StopIteration: # Time generator came to a halt
                        self.timegens.remove(tgen)
                        print('depleted')
                        if(len(self.timegens)==0): # Check if there are no more lapse to perform
                            loop = False
                            if self.call: # Alert the schedule to leap onto the next task
                                print('Calling back, last chance')
                                loop = self.call.alert()
                            if not loop : # No more events, thread is stopped
                                print('No more events')
                                return
                if loop: print('looping'); continue # Start over with the lapses
                dt = min(dts) 
                self.move(dt) # Advance with the smallest time lapse 
                self.clock += dt
                self.timeline.append(self.clock)
                if self.call : # Ask to check for new events
                    self.call.notify(dt)

        def move(self, dt):
            ''' Asks all animated objects to move with the tick dt'''
            if(dt!=0):
                for obj in self.animate:
                        obj.clocktick(dt)

        def affix(self, thread):
            ''' Join the thread = pause until the thread finishes

            Arguments
            ---------
            thread : Thread object - the thread that takes priority over this one 
            '''
            thread.join()

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
    
    def plan(self, eventname, func, *args, time = 0, duration = iter([0]), **kwargs):
        ''' Organise an event 

        Arguments
        ---------
        eventname : string - name displayed by the event once triggered
        func : function - function called by the event
        --------- multiple args
        *args : arguments passed onto the function of the event
        --------- Keyword args
        time : time at which the event must occur
        duration : iterable - indicates for how long must the event span
        **kwargs : passed onto the event
        '''
        evt = Event(eventname,duration,func,*args, **kwargs)
        if(time in self.events):
            self.events[time].append(evt)
        else:
            self.events[time] = [evt]
            self.times.append(time)
        self.times.sort()

    def urge(self, eventname, func, *args, delay = 0, duration = iter([0]), **kwargs):
        ''' Plan an event immediately, give an optional delay before it launches'''
        self.plan(eventname, func, *args, time = self.record+delay, duration = duration, **kwargs)

    def at(self, time):
        ''' Updates the schedule time and return the events due'''
        self.record = time
        if(len(self.times)!=0 and time>=self.times[0]):
            return self.events[self.times.pop(0)]
        else:
            return []

    def next(self):
        ''' Look for the next event due

        Returns
        -------
        the time that the next first event should take place
        '''
        if len(self.times)!=0: return self.times[0]
        return None

class Event:
    ''' Object containing a name, a task function and a slicing time generator'''

    def __init__(self, name, timegen, func, *args, **kwargs):
        self.name = name
        self.func = func
        self.args = args
        self.timegen = timegen
        self.type = kwargs.pop("type",False)
        self.key = kwargs.pop("key",self.name)

    def trigger(self):
        print(self.name)
        if(self.func is not None):
            self.func(*self.args)   