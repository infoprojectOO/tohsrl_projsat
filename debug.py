# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:59:17 2017

@author: ASiapan
"""
import six, subprocess, shlex, os, pathlib
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
import sys, string, os, subprocess
from scipy.optimize import curve_fit
#from satclass import *

#Input conversion
mpl.close("all")
fig, ax = mpl.subplots()
#ax = fig.add_subplot(111,projection = '3d')
# sat = ale_sat
# proj = projectile
#board = ax
#cam = ale_sat.camera
#atmosphere = earth.atm
#atm = earth.atm
# orbit = Orbit((a_sat,e_sat,i_sat,Om_sat,wp_sat))
def markplot(ax,x,y,id,color = 'r'):
    xi = [x[i] for i in id]
    yi = [y[i] for i in id]
    ax.scatter(xi,yi,c = color, marker = 'x', s = 50)
    
y = (np.linalg.norm(ghostectile.traj,axis=1)-earth.R)*0.001
x = ghostectile.time
marks = ghostbox.markindex

markplot(ax,x,y,marks)



    

