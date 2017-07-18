# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:26:28 2017

@author: ASiapan
"""

# DSMC code implementation

import math as m
import numpy as np
import random as rand
from itertools import product
from dsmcgraph import showfield

class Glocule:

    def __init__(self,r_gloc,v_gloc,cell):
        self.r = r_gloc
        self.v = v_gloc
        self.cell = cell

    def clocktick(self, dt):
        self.r += self.v*dt

class Cellule:

    def __init__(self, pos, size, index):
        self.pos = pos
        self.index = index
        self.size = size
        self.glocs = []

    def add_gloc(self, gloc):
        self.glocs.append(gloc)


n_molc = 1e+19 # /m³
MM_air = 0.0289644 # kg/mol
R_m = 8.315/MM_air # Jmol/kgK
U = 7000 # m/s
T = 200 # K

molc_scaling = 10**6
sigma_c = 1e-19 # m²
mfp = 1/(np.sqrt(2)*sigma_c*n_molc)

n_gloc = n_molc/molc_scaling

L_cell = mfp/100

x_bound = [-0.1,0.1]
y_bound = [-0.1,0.1]
z_bound = [-L_cell,L_cell]



def maxwell_distr(u_moy,T):
    s = np.sqrt(R_m * T)
    while(True):
        v_x = np.random.normal(u_moy[0],s)
        v_y = np.random.normal(u_moy[1],s)
        v_z = np.random.normal(u_moy[2],s)
        yield np.array((v_x,v_y,v_z))

def gridivide(xb,yb,zb,L):
    cells = []
    cell_indexes = []
    nx,ny,nz = m.ceil((xb[-1]-xb[0])/L),m.ceil((yb[-1]-yb[0])/L),m.ceil((zb[-1]-zb[0])/L)
    grid = (nx,ny,nz)

    L_size = ((xb[-1]-xb[0])/nx,(yb[-1]-yb[0])/ny,(zb[-1]-zb[0])/nz)
    xl = np.linspace(xb[0],xb[-1],nx+1)
    yl = np.linspace(yb[0],yb[-1],ny+1)
    zl = np.linspace(zb[0],zb[-1],nz+1)

    xc = np.linspace(xb[0]+L_size[0]*0.5,xb[-1]-L_size[0]*0.5,nx)
    yc = np.linspace(yb[0]+L_size[1]*0.5,yb[-1]-L_size[1]*0.5,ny)
    zc = np.linspace(zb[0]+L_size[2]*0.5,zb[-1]-L_size[2]*0.5,nz)
    cell_pos = np.array(list(product(xc,yc,zc)))
    cell_indexes = np.reshape(np.array([x for x in range(cell_pos.shape[0])]),(len(xc),len(yc),len(zc)))
    for c in range(cell_pos.shape[0]):
        cells.append(Cellule(np.array(cell_pos[c,:]),prod(L_size),c))
    np.reshape(cells,(len(xc),len(yc),len(zc)))
    return cells, cell_indexes, L_size, grid

cells, cell_indexes, cell_size, grid = gridivide(x_bound,y_bound,z_bound,L_cell)

def genmolc(cells,cell_size,n_gloc,vel_gen):
    glocs = []
    N_gloc_cell = int(np.round(prod(cell_size)*n_gloc))
    sample_glocs = []
    for i in range(N_gloc_cell):
        print('Gloc n° : ',i,'/', N_gloc_cell)
        r_gloc = np.array(cell_size)*np.array((rand.randrange(0,1)-0.5,rand.randrange(0,1)-0.5,rand.randrange(0,1)-0.5))
        v_gloc = next(vel_gen)
        sample_glocs.append((r_gloc,v_gloc))
    for c in np.ravel(cells):
        print('Cell n° : ',c)
        r_cell = c.pos
        for i in range(N_gloc_cell):
            r_gloc, v_gloc = sample_glocs[i]
            r_gloc = r_gloc + r_cell
            glocs.append(Glocule(r_gloc,v_gloc,c))
    return glocs

glocs = genmolc(cells,cell_size,n_gloc,maxwell_distr((U,0,0),T))

showfield(cells,(x_bound,y_bound,z_bound))