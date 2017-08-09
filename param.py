# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:07:51 2017

@author: ASiapan
"""

# Physical Parameters file

import math as m
from enum import Enum

# Earth

mu_G = 398600.4418 * 10**9 # m³/s²
R_G = 6378.1 * 10**3 # m 

# Atmosphere

MM_air = 0.0289644 # kg/mol
N_A = 6.023 * 10**23
R = 8.315
sigma_c = 1e-19 # m² mean collisionnal cross section
gamma = 1.4
b_air = 1.458e-6 # kg/m.s(K)^0.5
S_air = 110.4 # K

# Projectile parameters

d_p = 0.01 # m - diameter
S_p = d_p**2 * m.pi # m² - surface
# rho_p = 8.96 * 1000 # kg/m³ - material (copper)
# k_c = 401 # W/(m K)
# m_p = m.pi/6*rho_p*d_p**3 # kg - mass

# Physics

sigma_r = 5.670373 * 10**(-8) # W/m²K^4
k_B = 1.38064852*10**(-23) #J/K

# Satellite orbital parameters
a_sat = R_G + 350 * 10**3 # m - semimajor axis 
e_sat = 0.2 # ellipticity
i_sat = 0. * m.pi/180 # rad - inclination
Om_sat = 0. * m.pi/180 # rad - ascend node
wp_sat = 0. * m.pi/180 # rad - periapsis argument
nu0_sat = 0. * m.pi/180 # rad - starting anomaly


class Material(Enum):
	COPPER = (8.96 * 1000, 385,  401)
	ZINC = (8.96 * 1000, 385,  401)
	IRON = (8.96 * 1000, 385,  401)

	def __init__(self,rho, cal_cap, cal_cond):
		self.rho = rho
		self.C_cal = cal_cap # J/(kg·K)
		self.k_cal = cal_cond

	@property
	def mass_density(self):
		return self.rho

	@property
	def caloric_capacity(self):
		return self.C_cal

	@property
	def caloric_conductivity(self):
		return self.k_cal



