# Rarefied gas simulation coding Monte-Carlo

import math as m
import random as rand
from density import density
import matplotlib.pyplot as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def N_coll(n_mean):
	if(n_mean >= 50):
		n_coll = int(np.random.normal(n_mean,1/m.sqrt(n_mean)))
	else:
		n_coll = 0
		for i in range(int(n_mean)*4):
			n_coll += int(rand.random()>=0.5)
	return n_coll


def sphe2cart(r,lat,lon):
	x = r*np.cos(lat)*np.cos(lon)
	y = r*np.cos(lat)*np.sin(lon)
	z = r*np.sin(lat)
	return np.array([x,y,z])

# Environment parameters

mu_G = 398600 *1000**3 # m³/s²
R_G = 6378 * 1000 # m

alt_ini = 400 * 1000 # m
alt_tar = 100 * 1000 # m
MM_air = 0.0289644 # kg/mol
N_A = 6.023 * 10**23

userinput = False

# Location and time coordinates
if(userinput):
	latlong = input('ALE Latitude , Longitude : ').split(',')
	latitude = float(latlong[0])
	longitude = float(latlong[1])
else:
	longitude = 0 * m.pi/180 
	latitude = 90 * m.pi/180

a_dga = R_G+(alt_ini-alt_tar)*0.5 # half big axis
r_apog = alt_ini + R_G #starting

# Projectile parameters

d_p = 0.002 # m diameter
S_p = d_p* m.pi
rho_p = 8.96 * 1000 # kg/m³ copper
m_p = m.pi/6.*rho_p*d_p**3 # kg mass

Kn = 10


sigma_c = 1e-19 # m²

# Variables : initial conditions

#r_p = np.array([0.,0.,alt_ini+R_G]) # m
r_p = sphe2cart(r_apog,latitude,longitude) # m
#r_p = r_sat
v_p = np.array([7000.,0.,0.]) # m/s
#v_p = 0.9*v_sat
#v_p[0] = m.sqrt(2*mu_G*((1/r_apog)-1/(2*a_dga)))
#v_p[0] = m.sqrt(2*mu_G*(1/(2*r_apog))) # Circular Orbit at 400 km
a_p = np.array([0.,0.,-mu_G/(np.linalg.norm(r_p))**2])

N_itmax = 10**4

dt = 10**(-3)*R_G/np.linalg.norm(v_p)

# Regime
rho = density(100)

output = atm.nrlmsise_output()
inputatm = atm.nrlmsise_input(year = 2000, doy = 1, sec = 0.0, alt=100.0, g_lat=30., g_long=0.0,
                 lst=0.0, f107A=150., f107=150.0, ap=4.0, ap_a=None)
flags = atm.nrlmsise_flags()
flags.switches[0] = 1
gtd7(inputatm,flags,output)
rho2 = output.d[5]
T = output.t[1]
b_air = 1.458e-6 # kg/ms(K)^0.5
S_air = 110.4 # K
visc = b_air*T**(1.5)/(T+S_air) # Sutherland for air


Re = rho2*np.linalg.norm(v_p)**2*d_p/visc

# Storage

t = [0]
speed_p = v_p[:]
traj_p = r_p[:]
n_it = 0



while (Kn >= 1 and n_it <= N_itmax):

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
 
	volume = S_p*np.linalg.norm(v_p)*dt # Tube volume covered by projectile during dt
	n_volmean = density(alt*0.001)*N_A/(MM_air)
	lpm = 1/(m.sqrt(2)*n_volmean*sigma_c)
	n_mean = n_volmean*volume
	n_coll = N_coll(n_mean)
	dv_coll = - n_coll*MM_air/(N_A*m_p) * v_p
	
	a_p = -r_p/np.linalg.norm(r_p)*mu_G/(np.linalg.norm(r_p))**2
	v_p_old = v_p
	v_p = v_p + a_p*dt + dv_coll

	#r_p = r_p + (v_p_old+v_p)*0.5*dt
	r_p = r_p + (2*v_p)*0.5*dt
   alt = np.linalg.norm(r_p)-R_G


	traj_p = np.vstack((traj_p,r_p))
	speed_p = np.vstack((speed_p,v_p))
	t = np.vstack((t,t[-1]+dt))

	Kn = lpm/d_p
	n_it += 1


# Plot data

new_Y = traj_p[0,:]/np.linalg.norm(traj_p[0,:])
new_X = speed_p[0,:]/np.linalg.norm(speed_p[0,:])
traj_projection = traj_p.dot(np.vstack((new_X,new_Y)).transpose())


# Trajectory in osculatory plane 
fig = mpl.figure(1)
# earth = mpl.pyplot.Circle(0,0,R_G,color = 'b')
# ax = fig.gca()
# ax.add_artist(earth)
mpl.plot(R_G*np.sin(np.linspace(0,2*m.pi,num=1000)),R_G*np.cos(np.linspace(0,2*m.pi,num=1000)),'b-')
mpl.plot(traj_projection[:,0],traj_projection[:,1],'r-')

# Altitude decrease with time 
mpl.figure(2)
mpl.plot(t,(np.linalg.norm(traj_p,axis=1)-R_G)*0.001,'r-')
mpl.xlabel('Time (s)')
mpl.ylabel('Altitude (km)')

mpl.figure(3)
mpl.plot((np.linalg.norm(traj_p,axis=1)-R_G)*0.001,np.linalg.norm(speed_p,axis=1)*0.001,'r-')
mpl.xlabel('Altitude (km)')
mpl.ylabel('Speed (km/s)')

# Earth plotting sphere
fig = mpl.figure(4)
ax = fig.add_subplot(111, projection='3d')
phi = np.linspace(0, 2 * np.pi, 10)
theta = np.linspace(0, np.pi, 10)
xg = R_G * np.outer(np.cos(phi), np.sin(theta))
yg = R_G * np.outer(np.sin(phi), np.sin(theta))
zg = R_G * np.outer(np.ones(np.size(phi)), np.cos(theta))
ax.plot_wireframe(xg, yg, zg)
ax.plot(traj_p[:,0], traj_p[:,1], traj_p[:,2],color = 'r')


