# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:59:17 2017

@author: ASiapan
"""
import six, subprocess, shlex, os, pathlib
from enum import Enum
import atmosphere.nrlmsise_00_header as atm
from atmosphere.nrlmsise_00 import gtd7
from six.moves import map, xrange, zip, reduce
import numpy as np
import matplotlib.pyplot as mpl
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import colors as mcolors
import methutil as methu
import sys, string, os, subprocess
from scipy.optimize import curve_fit
import dsmc
import glob
import graphics as gp
from olog import olog
from cycler import cycler
#from satclass import *

#Input conversion
# mpl.close("all")
# fig, ax = mpl.subplots()
#ax = fig.add_subplot(111,projection = '3d')
# sat = ale_sat
# proj = projectile
#board = ax
#cam = ale_sat.camera
#atmosphere = earth.atm
#atm = earth.atm
# # orbit = Orbit((a_sat,e_sat,i_sat,Om_sat,wp_sat))
# hypbox = mechanics.boxes[projectile]
# ghostectile = hypbox.ghost
# ghostbox = hypbox.ghostbox
# dsim = hypbox.dsmc
sparta_folder = "\\dsmc\\"
name = "simi"

aeroreg = lastghostreg

Ma2 = aeroreg.Ma**2
gamma = aeroreg.air.gamma
T_s = (2*gamma*Ma2-(gamma-1))*((gamma-1)*Ma2+2)/((gamma+1)**2*Ma2)*aeroreg.T
if(T_s)>8000: T_s = 8000
rho_s = aeroreg.air.rho*(gamma+1)*Ma2/((gamma-1)*Ma2 + 2)
S = ((gamma-1)*Ma2 + 2)/(2*gamma*Ma2-(gamma-1))*np.sqrt(aeroreg.T*gamma/(T_s*2))
v_therm = aeroreg.air.v_therm*np.sqrt(T_s/aeroreg.T)*np.sqrt(m.pi)/2
N_i = aeroreg.S*(rho_s/aeroreg.air.MM_air*aeroreg.air.N_A)*v_therm*((S**2+0.5)*m.erf(S)/S + m.exp(-S**2)/m.sqrt(m.pi))
