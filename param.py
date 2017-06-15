# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:07:51 2017

@author: ASiapan
"""

# Physical Parameters file

import math as m

# Earth

mu_G = 398600 * 10**9 # m³/s²
R_G = 6378 * 10**3 # m 

# Atmosphere

# Projectile parameters

d_p = 0.002 # m - diameter
S_p = d_p* m.pi # m² - surface
rho_p = 8.96 * 1000 # kg/m³ - material (copper)
m_p = m.pi/3*rho_p*d_p**3 # kg - mass


# Satellite orbital parameters
a_sat = R_G + 400 * 10**3 # m - semimajor axis 
e_sat = 0.2 # ellipticity
i_sat = 0. * m.pi/180 # rad - inclination
Om_sat = 0. * m.pi/180 # rad - ascend node
wp_sat = 0. * m.pi/180 # rad - periapsis argument
nu0_sat = 0. * m.pi/180 # rad - starting anomaly