# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:07:51 2017

@author: ASiapan

This file regroups all the numerical parameter values used in the code.
They are grouped together according to the objects they apply on
"""

# Physical Parameters file

import math as m
from enum import Enum
import datetime as dt

# Physics

sigma_r = 5.670373 * 10**(-8) # W/m²K^4
k_B = 1.38064852*10**(-23) #J/K

# Earth

mu_G = 398600.4418 * 10**9 # m³/s² - Earth gravity constant
R_G = 6378.137 * 10**3 # m - radius at the equator
R_pole = 6356.752 * 10**3 # m - radius to the poles
sidereal_G = 86164.1 # s - sidereal day duration
solarday_G = 86400 # s- solar day duration
tilt_G = 23.44*m.pi/180 # rad - earth tilt angle

# Atmosphere

MM_air = 0.0289644 # kg/mol
N_A = 6.023 * 10**23
R = k_B*N_A # K/(J.mol)
R_m = R/MM_air #K/(kg.J)
sigma_c = 1e-19 # m² - mean collisionnal cross section
gamma = 1.4
b_air = 1.458e-6 # kg/m.s(K)^0.5
S_air = 110.4 # K

# Projectile parameters

d_p = 0.01 # m - diameter
S_p = d_p**2 * m.pi # m² - surface
T_p = 200 # K - projectile temperature

# Satellite orbital parameters
a_sat = R_G + 350 * 10**3 # m - semimajor axis 
e_sat = 0.2 # ellipticity
i_sat = 0. * m.pi/180 # rad - inclination
Om_sat = 0. * m.pi/180 # rad - ascend node
wp_sat = 0. * m.pi/180 # rad - periapsis argument
nu0_sat = 0. * m.pi/180 # rad - starting anomaly


class Material(Enum): # atomic mass, density, caloric capacity, caloric conductivity, melting temperature, fusion specific heat
    COPPER = (63.546, 8.96 * 1000, 385,  401, 1357.77, 13260/0.063546, 'Cu')
    ZINC = (65.38, 7.134 * 1000, 388,  116, 692.68, 7320/0.06538, 'Zn')
    IRON = (55.845, 7.874 * 1000, 449, 80.4, 1811, 13810/0.055845, 'Fe')
    TITANIUM = (47.867, 4.506 * 1000, 524, 21.9, 1941, 14150/0.047867, 'Ti')

    def __init__(self,ma,rho, cal_cap, cal_cond, T_melt, h_melt, id = ''):
        self.ma = 0.001*ma # kg/mol
        self.rho = rho # kg/m³
        self.C_cal = cal_cap # J/(kg·K)
        self.k_cal = cal_cond # W/(m K)
        self.T_melt = T_melt # K
        self.h_melt = h_melt # J/kg
        self._id = id

    @property
    def mass_density(self):
        return self.rho

    @property
    def caloric_capacity(self):
        return self.C_cal

    @property
    def caloric_conductivity(self):
        return self.k_cal

    @property
    def mol_mass(self):
        return self.ma

    @property
    def melting_T(self):
        return self.T_melt

    @property
    def melt_specheat(self):
        return self.h_melt

    @property
    def id(self):
        return self._id


class Reference:
    ''' Convenient class for regrouping time and frame references of the code '''
    _refdate = dt.datetime(2000,1,1,12)

    @classmethod
    def get_refdate(cls):
        ''' Returns the date taken as the default time-frame reference'''
        return Reference._refdate

    @classmethod
    def get_refloc(cls):
        ''' Returns the geographic coordinates of a place of reference on the globe'''
        return 35.6895*m.pi/180, 139.69171*m.pi/180, 'Tokyo' # Tokyo coordinates : Latitude, Longitude

