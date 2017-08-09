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
import os
from enum import Enum
import pathlib
import glob

sparta_folder = "\\dsmc\\"

class SimInput(Enum):
    AIR_NDENSITY = "nv_air"
    AIR_TEMPERATURE = "T_air"
    SPEED = "v_air"
    PROJ_TEMPERATURE = "T_proj"

class Sparta_Writer:
    base_name = "axi.in"
    example_input = "axi.ex"
    prime_input = "axi.init"
    input_path = os.getcwd()+sparta_folder
    token = "§"

    def __init__(self):
        self.clear_folder()

    def write_simulation(self,sim_param,fresh = False):
        simid = sim_param.pop("name")
        simname, siminstr = self.prepare_simfile(simid,fresh)
        with open(self.input_path+simname,"w") as f:
            for k,v in sim_param.items():
                siminstr = siminstr.replace(self.token+k.value,str(v),1)
            f.write(siminstr)
        sim_param["simname"] = simname
        return simname


    def prepare_simfile(self,id,fresh):
        path = self.input_path
        protofile = self.example_input
        filename = self.base_name + "_" + id
        if fresh: 
            protofile = self.prime_input
            filename = self.prime_input + "_" + id
        with open(path+protofile,"r", encoding='utf-8') as sample_file:
            sample_instr = sample_file.read()
        simfile = pathlib.Path(path+filename)
        simfile.write_text("")
        return filename, sample_instr

    def clear_folder(self):
        clearlist = glob.glob("."+sparta_folder+self.base_name+"_*")
        clearlist.extend(glob.glob("."+sparta_folder+self.prime_input+"_*"))
        print("Clearing folder for new simulation inputs : ",str(len(clearlist))," file(s) removed")
        for clearfile in clearlist:
            os.remove(clearfile)






'''
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

molc_scaling = 10**9
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
        cells.append(Cellule(np.array(cell_pos[c,:]),np.prod(L_size),c))
    np.reshape(cells,(len(xc),len(yc),len(zc)))
    return cells, cell_indexes, L_size, grid

cells, cell_indexes, cell_size, grid = gridivide(x_bound,y_bound,z_bound,L_cell)

def genmolc(cells,cell_size,n_gloc,vel_gen):
    glocs = []
    N_gloc_cell = int(np.round(np.prod(cell_size)*n_gloc))
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

# glocs = genmolc(cells,cell_size,n_gloc,maxwell_distr((U,0,0),T))

# showfield(cells,(x_bound,y_bound,z_bound))
'''