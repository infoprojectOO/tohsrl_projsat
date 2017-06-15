# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 18:35:46 2017

@author: ASiapan
"""

import numpy as np
import matplotlib.pyplot as mpl
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation



def draw1(num):
    # NOTE: there is no .set_data() for 3 dim data...
    line1.set_data(r_proj[0,0:num],r_proj[1, 0:num])
    line1.set_3d_properties(r_proj[2, 0:num])
    return line1

def init():
    return [line1]

# Attaching 3D axis to the figure
fig = mpl.figure()
ax = p3.Axes3D(fig)

# Fifty lines of random 3-D lines
#data = 0.04*np.array(range(0,25,1))*np.ones((3,25))
r_proj = 0.04*np.array([list(range(25)),list(range(25)),list(range(25))])

# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
line1, = ax.plot([], [], [])

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, draw1, 25, interval=50, blit=True, init_func = init)

plt.show()