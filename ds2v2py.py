# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 12:23:47 2017

@author: ninja_000
"""
import math as m
import numpy as np
from satclass import *
import random as rand

# Environment parameters

mu_G = 398600 *1000**3 # m³/s²
R_G = 6378 * 1000 # m

alt_ini = 400 * 1000 # m
alt_tar = 100 * 1000 # m
MM_air = 0.0289644 # kg/mol
N_A = 6.023 * 10**23


def N_coll(n_mean):
	if(n_mean >= 50):
		n_coll = int(np.random.normal(n_mean,1/m.sqrt(n_mean)))
	else:
		n_coll = 0
		for i in range(int(n_mean)*4):
			n_coll += int(rand.random()>=0.5)
	return n_coll

def manSimIter(proj, earth,dragc):
    
#    Constant
    gamma = 1.4 #assume perfect air
    n_it = 0
    N_itmax = 10**4
#    sigma_c = 1e-19 # m²
    t= [0]
    R = 286 #m^2 /s^2/K
   
    v_p = proj.v
    r_p = proj.r
    m_p = proj.m
    S_p = proj.S
    lat, lon = earth.pos2coord(r_p)
    air = earth.atm.at(0,np.linalg.norm(r_p),lat, lon)
    
    dt = 10**(-3)*R_G/np.linalg.norm(v_p)
    speed_p = v_p[:]
    traj_p = r_p[:]
    a_v = m.sqrt(gamma*R*air.T)
    mach = np.linalg.norm(v_p)/a_v
    alt_ini = np.linalg.norm(r_p)-R_G
    alt = np.linalg.norm(r_p)-R_G
#    print(mach)
   
    while (alt >= alt_ini-5000 and n_it <= N_itmax):
#         volume = S_p*np.linalg.norm(v_p)*dt # Tube volume covered by projectile during dt
#         n_volmean = density(alt*0.001)*N_A/(MM_air)
#         lpm = 1/(m.sqrt(2)*n_volmean*sigma_c)
#         n_mean = n_volmean*volume
#         n_coll = N_coll(n_mean)
#         dv_coll = - n_coll*MM_air/(N_A*m_p) * v_p
         lat, lon = earth.pos2coord(r_p)
         air = earth.atm.at(0,np.linalg.norm(r_p)-R_G,lat, lon)
         T = air.T
         rho = density(alt/1000)
         nv = air.nv
         dv_drag = - dragc*S_p*np.linalg.norm(v_p)*v_p*rho*dt/m_p
         a_p = -r_p/np.linalg.norm(r_p)*mu_G/(np.linalg.norm(r_p))**2
#         v_p_old = v_p
         a_v = m.sqrt(gamma*R*T)
         v_p += a_p*dt + dv_drag
         v_pnorm = np.linalg.norm(v_p)
         mach = v_pnorm/a_v
         #r_p = r_p + (v_p_old+v_p)*0.5*dt
         r_p = r_p + (2*v_p)*0.5*dt
         alt = np.linalg.norm(r_p)-R_G
         traj_p = np.vstack((traj_p,r_p))
         speed_p = np.vstack((speed_p,v_p))
         t = np.vstack((t,t[-1]+dt))
#         Kn = lpm/d_p
         n_it += 1
         g_a = mu_G/np.linalg.norm(r_p)**2
#         print(str(np.linalg.norm(dv_drag))+ ',' + str(g_a))
         
#         print('Iteration: ' + str(n_it))
#         print('at Alt = ' + str(alt))
    proj.v = v_p
    proj.r = r_p
    proj.m = m_p
    proj.S = S_p
#    print('Iteration: ' + str(n_it))
#    print('at Alt = ' + str(alt))
#    print('Mac = ' + str(mach))
    return rho, nv, T, v_pnorm, alt, proj
	
# alt = altitude (m)
	# r_p = projectile vector position (m)
	# R_G = earth radius (m)
	# S_p = projectile surface (m²)
	# N_A = avogadro number (for moles)
	# MM_air = air mean mass per mole (kg/mol)
	# sigma_c = air mean collisional cross section (m²)
	# m_p = projectile mass (kg)
	# v_p = projectile velocity vector (m/s)
	# mu_G = earth's gravity constant (m³/s²)
	# d_p = projectile diameter (m)
	# a_p = projectile gravity acceleration vector (m/s²)