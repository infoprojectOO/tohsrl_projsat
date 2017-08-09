# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 18:21:03 2017

@author: ASiapan
"""

# DSMC flowfield display

import matplotlib.pyplot as mpl
import numpy as np
from itertools import product

def showfield(cells,boundaries):
    matim = formMatrix(cells)
    fig  = mpl.figure()
    mpl.imshow(matim, extent=(min(boundaries[0]), max(boundaries[0]), min(boundaries[1]), max(boundaries[1])),interpolation='nearest',cmap=cm.gist_rainbow)
    points = np.array(list(product(range(3),range(4))))

    mpl.grid(True,xdata = np.linspace(-1,1,11), ydata = np.linspace(-1,1,11))
    mpl.show()


def formMatrix(cells):
    mat = np.zeros(cells.shape[0:2])
    for row in range(cells.shape[0]):
        for col in range(cells.shape[1]):
            density = 0
            for c in cells[row,col,:]:
                density += len(c.glocs)/c.size
            mat[col,row] = density

def showgrid(ax,grid,boundaries):
    x = np.linspace(*boundaries[0],grid[0]+1)
    y = np.linspace(*boundaries[1],grid[1]+1)

    for i in range(len(y)):
        ax.plot(x,[y[i]]*len(x),'r')
    for j in range(len(y)):
        ax.plot([x[j]]*len(y),y,'r')

