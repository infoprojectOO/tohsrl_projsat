# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:59:17 2017

@author: ASiapan
"""
import six
from density import density
import atmosphere.nrlmsise_00_header as atm
from atmosphere.nrlmsise_00 import gtd7
from six.moves import map, xrange, zip, reduce
import numpy as np
from matplotlib.axes import Axes, rcParams
import matplotlib.pyplot as mpl
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import colors as mcolors
from methutil import cube

#Input conversionµ
# mpl.close("all")
# fig = mpl.figure()
# ax = fig.add_subplot(111,projection = '3d')
# #ax = p3.Axes3D(fig)
# sat = ale_sat
# proj = projectile
# board = ax
# cam = ale_sat.camera
# rim = 0.5

r_sat = [0,0,0]

rho = density(100)

output = atm.nrlmsise_output()
inputatm = atm.nrlmsise_input(year = 2000, doy = 1, sec = 0.0, alt=100.0, g_lat=30., g_long=50.0,
                 lst=0.0, f107A=150., f107=150.0, ap=4.0, ap_a=None)
flags = atm.nrlmsise_flags()
flags.switches[1:] = [1]*(len(flags.switches)-1)

flags.switches[0] = 1
gtd7(inputatm,flags,output)
rho2 = output.d[5]
T = output.t[1]
b_air = 1.458e-6 # kg/ms(K)^0.5
S_air = 110.4 # K
visc = b_air*T**(1.5)/(T+S_air)
Re = rho2*7000**2*0.002/visc
print('First : ',rho2,'\n')

for i in range(1,len(flags.switches)):
    try:
        flags.switches[i] = 0
        flags.switches[i-1] = 1
        gtd7(inputatm,flags,output)
    except(ZeroDivisionError):
        print('error :', i)
    print('test :', i, (rho2-output.d[5])/rho2)