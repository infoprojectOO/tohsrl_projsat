# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:59:32 2017

@author: ASiapan
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')#(xlim=([-1.1*R_G, 1.1*R_G]), ylim=([-1.1*R_G, 1.1*R_G]), zlim = [-R_G*1.1,R_G*1.1])
ax.set_xlim([-1.1*R_G, 1.1*R_G])
ax.set_ylim([-1.1*R_G, 1.1*R_G])
ax.set_zlim([-1.1*R_G, 1.1*R_G])

line, = ax.plot([], [], [], lw=2)

r_proj = projectile.traj

# initialization function: plot the background of each frame
def init():
    return [line]

# animation function.  This is called sequentially
def draw(index):
    # update the data
    # t, r_sat, att_sat, r_proj = data
    line.set_data(r_proj[0:index,0],r_proj[0:index,1])
    line.set_3d_properties(r_proj[0:index,2])
    return [line]

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, draw, len(r_proj), init_func=init,
                                interval=20, blit=True, repeat = False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()