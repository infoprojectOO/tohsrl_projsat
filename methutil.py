# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:32:55 2017

@author: ASiapan
"""

#Tricky method file

import numpy as np
import math as m
from itertools import product
from mathutils import Quaternion, Vector
from density import density
from sympy.solvers.solveset import nonlinsolve as nls
from sympy.core.symbol import symbols
from sympy import exp, sqrt
import scipy.optimize as opt

def sphe2cart(r,lat,lon):
    x = r*np.cos(lat)*np.cos(lon)
    y = r*np.cos(lat)*np.sin(lon)
    z = r*np.sin(lat)
    return np.array([x,y,z])

def permute(array, line1, line2, axis):
    if(axis==0):
        temp = array[:,line1].copy()
        array[:,line1] = array[:,line2]
        array[:,line2] = temp
    elif(axis==1):
        temp = array[line1,:]
        array[line1,:] = array[line2,:]
        array[line2,:] = temp

def cube(position = np.zeros((3)), attitude = Quaternion(), size = 1):
    # combi_rel = np.array(list(product(r,r,r))).T
 #    combi_abs = m_rot.dot(combi_rel)
 #    bcombi = m_rot.dot(np.array(list(product(r,r,[0,0]))).T)[0:3,0:8:2]
 #    combi1,combi2 = np.hsplit(combi_abs,2)
 #    rows,cols = combi_abs.shape
 #    combis = np.array(list(zip(combi1,combi2))).reshape(rows*2,int(cols/2))
 #    xb = np.ones((2,2))*r_sat[0]+bcombi[0,:].reshape(2,2)
 #    yb = np.ones((2,2))*r_sat[1]+bcombi[1,:].reshape(2,2)
 #    zb = np.ones((2,2))*r_sat[2]+bcombi[2,:].reshape(2,2)
 #    x_sat = np.ones(combis.shape)*r_sat[0]+combis
 #    y_sat = np.ones(combis.shape)*r_sat[1]+np.roll(combis,-2,axis=0)
 #    z_sat = np.ones(combis.shape)*r_sat[2]+np.roll(combis,-4,axis=0)
    r = [-size*0.5,size*0.5]
    vortex = np.flipud(np.array(list(product(r,r,r))).T)
    permute(vortex,2,3,axis=0)
    permute(vortex,6,7,axis=0)
    vxy1,vxy2 = np.hsplit(vortex,2)
    vzx1,vzx2 = np.hsplit(np.roll(vortex,-1,axis=0),2)
    vyz1,vyz2 = np.hsplit(np.roll(vortex,-2,axis=0),2)
    #vxy1 = np.hstack([vxy1,vxy1[:,0]])
    vyz2 = np.roll(vyz2,-1,axis=1)
    vzx1 = np.fliplr(vzx1)
    cubepts = np.hstack([vxy1,vzx2,vyz2,vzx1,vyz1,vxy2])

    m_rot = np.array(attitude.to_matrix())
    cubeabs = m_rot.dot(cubepts) + position.reshape((3,1))
    return cubeabs

def plane_proj(vec1,vec2):
    vec1 = vec1/np.linalg.norm(vec1)
    vec2 = vec2-vec1.dot(vec2)*vec1
    vec2 = vec2/np.linalg.norm(vec2)
    vec3 = np.cross(vec1,vec2)
    return np.vstack((vec1,vec2,vec3))

def angle_vec(vec1,vec2,ref):
    nan = float('nan')
    # ax1 = int(len(vec1.shape)>=2)
    # ax2 = int(len(vec2.shape)>=2)
    if(vec1.ndim==1):
        vec1 = vec1.reshape((1,vec1.size))
        vec2 = vec2.reshape((1,vec1.size))
    n = vec1.shape[0]
    vec1 = vec1/np.linalg.norm(vec1,axis=1).reshape((n,1))
    vec2 = vec2/np.linalg.norm(vec2,axis=1).reshape((n,1))
    print(vec1[vec1==nan],vec2[vec2==nan])
    vec1[vec1==nan] = 0 # Zero division fixing
    vec2[vec2==nan] = 0 # Zero division fixing

    vec3 = np.cross(vec1,vec2)
    sign = np.sign(np.sum(vec3*ref,axis = 1))
    sin = sign*np.linalg.norm(vec3,axis = 1)
    cos = np.sum(vec1*vec2,axis = 1)
    cos[cos>1.] = 1.
    cos[cos<1.] = 1.
    ang = angs = np.arcsin(sin)
    angc = np.arccos(cos)
    for i in range(len(ang)):
        if(angs[i]<=0 and angc[i]>=np.pi*0.5):
            ang[i] = 2*m.pi-angc[i]
        elif(angs[i]>=0 and angc[i]>=np.pi*0.5):
            ang[i] = angc[i]
    # elif(angs[i]<=0 and angc[i]<=np.pi*0.5):
    #   ang = angs
    return ang#*180/np.pi

def plane_intersection(n1,n2,diref=np.array([1,0,0])):
    n1 = n1/np.linalg.norm(n1)
    n2 = n2/np.linalg.norm(n2)
    n3 = np.cross(n1,n2)
    if(sum(n3)<=1e-16):
        n3 = diref
    else:
        n3 = n3/np.linalg.norm(n3)
    return n3

def fit_atm(atmosphere):
    lat, lon, name = Reference.get_refloc()
    rho_scale, h_scale = atmosphere.profile(np.linspace(50,100,801),lat, lon)
    m1 = int(len(rho_scale)*1/4)
    m2 = int(len(rho_scale)*3/4)
    slope = (np.log(rho_scale[m2])-np.log(rho_scale[m1]))/(h_scale[m2]-h_scale[m1])
    def expo(x,const,slope):
        return const*np.exp(slope*x) # + offset not taken into account
    params, s = opt.curve_fit(expo, h_scale,rho_scale, (rho_scale[0],slope))
    # mpl.figure()
    # mpl.plot(h_scale,rho_scale,'b-')
    # mpl.plot(h_scale,expo(h_scale,*params))
    return params

def integrate(y,x):
    dx = np.roll(x,-1,axis=0)-x
    if(dx.ndim>1):
        dx = np.linalg.norm(dx,axis=1)
    dx = dx[:-1]
    inty = []
    for yp in y:
        intup = np.vdot(yp[:-1],dx) # Prograde rectangle integration
        intdo = np.vdot(yp[1:],dx) # Retrograde rectangle integration
        inty.append((intup+intdo)/2)
    return inty

def minimum(dic):
    dico = dic.copy()
    km, vm = dico.popitem()
    dico[km] = vm
    for k,v in dico.items():
        if(v<vm):
            km = k
            vm = v
    return (km)

def minima(dico,n):
    if(len(dico.values())<n):
        raise Exception("More indexes required than number of elements")
    elif(n<=0):
        raise Exception("Unacceptable number of minima search")
    if n==1:
        return [minimum(dico)]
    else:
        prec = minimum(dico)
        copy = dico.copy()
        copy[prec] = max(dico.values())+1
        return [prec]+minima(copy,n-1)






    