# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 10:59:16 2017

@author: ASiapan
"""
import numpy as np
import os
import PIL
from satclass import *
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import axes3d, proj3d
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3D
from matplotlib.patches import FancyArrowPatch
from mathutils import Quaternion
import methutil as methu
from methutil import cube

R_G = 6378 * 10**3 # m 
o = Orbit((30,0,0,0,0))

def convert2vert(X,Y,Z,*args,**kwargs):
    Z = np.atleast_2d(Z)
    # TODO: Support masked arrays
    X, Y, Z = np.broadcast_arrays(X, Y, Z)
    rows, cols = Z.shape

    rstride = 1
    cstride = 1

    polys = []

    for rs in range(0, rows-1, rstride):
        for cs in range(0, cols-1, cstride):
            ps = []
            for a in (X, Y, Z) :
                ztop = a[rs,cs:min(cols, cs+cstride+1)]
                zleft = a[rs+1:min(rows, rs+rstride+1),
                          min(cols-1, cs+cstride)]
                zbase = a[min(rows-1, rs+rstride), cs:min(cols, cs+cstride+1):][::-1]
                zright = a[rs:min(rows-1, rs+rstride):, cs][::-1]
                z = np.concatenate((ztop, zleft, zbase, zright))
                ps.append(z)

            # The construction leaves the array with duplicate points, which
            # are removed here.
            ps = list(zip(*ps))
            lastp = np.array([])
            ps2 = [ps[0]] + [ps[i] for i in range(1, len(ps)) if ps[i] != ps[i-1]]
            avgzsum = sum(p[2] for p in ps2)
            polys.append(ps2)
    return polys

class satPainter:
    # Painter.add(Satellite,satPainter)
    
    def __init__(self,sat):
        self.sat = sat
        self.rim = sat.width
        self.r = sat.orbit.getPos(sat.nu)
        self.att = sat.attitude
        self._cube = cube(np.zeros((3)),Quaternion(),self.rim)
        self._cam = self.camInit(self.sat.camera.width,-0.5*self.sat.width)
        ptscube, ptscam = self.center(self.r,self.att)
        xs,ys,zs = ptscube[0,:].reshape(12,2), ptscube[1,:].reshape(12,2), ptscube[2,:].reshape(12,2)
        xc,yc,zc = ptscam[0,:].reshape((26,26)), ptscam[1,:].reshape((26,26)), ptscam[2,:].reshape((26,26)) 
        verts = convert2vert(xs,ys,zs)
        self.body = Poly3DCollection(verts, color = 'm')
        vertscam = convert2vert(xc,yc,zc)
        self.cam = Poly3DCollection(vertscam, color = 'c')

    def camInit(self,width,z_offset):
        xcam,ycam = np.meshgrid(np.linspace(-0.5*width,0.5*width,26),np.linspace(-0.5*width,0.5*width,26))
        print(z_offset)
        zcam = -np.sqrt(xcam**2+ycam**2)+z_offset
        return np.vstack((xcam.reshape((1,xcam.size)),ycam.reshape((1,ycam.size)),zcam.reshape((1,zcam.size))))
    
    def center(self,pos,att):
        m_rot = np.array(att.to_matrix())
        cubeabs = m_rot.dot(self._cube) + pos.reshape((3,1))
        camabs = m_rot.dot(self._cam) + pos.reshape((3,1))
        return cubeabs, camabs
        
    def update(self, index):
        self.index = index
        self.r = self.sat.traj[index,:]
        self.v = self.sat.speed[index,:]
        self.att = self.sat.attitudes[index]
        pts, ptscam = self.center(self.r,self.att)
        xs,ys,zs = pts[0,:].reshape(12,2), pts[1,:].reshape(12,2), pts[2,:].reshape(12,2)
        xc,yc,zc = ptscam[0,:].reshape((26,26)), ptscam[1,:].reshape((26,26)), ptscam[2,:].reshape((26,26))
        self.body.set_verts(convert2vert(xs,ys,zs))
        self.cam.set_verts(convert2vert(xc,yc,zc))
        return self.get_artists()

    def get_artists(self):
        return [self.body, self.cam]

        # Plot satellite velocity at final position
        # mag = 200
        # xs,ys,zs = r_sat
        # xv,yv,zv = mag*v_sat+r_sat
        # varrow = Arrow3D([xs,xv],[ys,yv],[zs,zv], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
        # board.add_artist(varrow)

        # # Plot satellite Axis
        # xd, yd, zd = m_rot.dot(np.array([0,0,sat.width*2]))+r_sat
        # sarrow = Arrow3D([xs,xd],[ys,yd],[zs,zd],mutation_scale=20, lw=1, arrowstyle="-|>", color="y")
        # board.add_artist(sarrow)

class earthPainter:
    # Painter.add(Earth,earthPainter)

    def __init__(self,earth, globe = False):
        self.earth = earth
        phi = np.linspace(0, 2 * np.pi, 10)
        theta = np.linspace(0, np.pi, 10)
        xg = R_G * np.outer(np.cos(phi), np.sin(theta))
        yg = R_G * np.outer(np.sin(phi), np.sin(theta))
        zg = R_G * np.outer(np.ones(np.size(phi)), np.cos(theta))
        ax = mpl.figure(frameon = False).add_subplot(111,projection = '3d')
        self.frame = ax.plot_wireframe(xg,yg,zg)
        self.frame.remove()
        ax.get_figure().clear()
        

        self.phi0 = earth.rotang
        self.rotation = Quaternion([0,0,1],self.phi0)
        self.xm = R_G * np.cos(self.phi0)*np.sin(theta)
        self.ym = R_G * np.sin(self.phi0)*np.sin(theta)
        self.zm = R_G * np.cos(theta)

        self.meridian = Line3D(self.xm,self.ym,self.zm,color = 'y')
        self._see_globe = globe
        if(globe):
            self.load_globe()


    def update(self,index):
        phirot = self.earth.phi[index]
        self.rotation.angle = phirot
        
        if(self._see_globe):
            m_rot = np.array(self.rotation.to_matrix())
            tmp = np.vstack((self.xglobe.reshape((1,self.xglobe.size)), self.yglobe.reshape((1,self.yglobe.size)), self.zglobe.reshape((1,self.zglobe.size))))
            xr,yr,zr = m_rot.dot(tmp)
            xr,yr,zr = xr.reshape(self.xglobe.shape), yr.reshape(self.xglobe.shape), zr.reshape(self.xglobe.shape)
            self.globe.set_verts(convert2vert(xr,yr,zr))

        m_rot = np.array(Quaternion([0,0,1],phirot-self.phi0).to_matrix())
        xl,yl,zl = m_rot.dot(np.vstack((self.xm,self.ym,self.zm)))
        self.meridian.set_data(xl,yl)
        self.meridian.set_3d_properties(zl)

        return self.get_artists()

    def get_artists(self):
        a = [self.meridian]
        if self._see_globe:
            a.append(self.globe)
        return a


    def load_globe(self):
         # load bluemarble with PIL
        bm = PIL.Image.open('cylmapearth.jpg')
        # it's big, so I'll rescale it, convert to array, and divide by 256 to get RGB values that matplotlib accept 
        bm = np.array(bm)/256.

        # coordinates of the image 
        lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
        lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 

        self.xglobe = R_G*np.outer(np.cos(lons), np.cos(lats)).T
        self.yglobe = R_G*np.outer(np.sin(lons), np.cos(lats)).T
        self.zglobe = R_G*np.outer(np.ones(np.size(lons)), np.sin(lats)).T
        self.colormap = bm

        
        m_rot = np.array(self.rotation.to_matrix())
        tmp = np.vstack((self.xglobe.reshape((1,self.xglobe.size)), self.yglobe.reshape((1,self.yglobe.size)), self.zglobe.reshape((1,self.zglobe.size))))
        xr,yr,zr = m_rot.dot(tmp)
        xr,yr,zr = xr.reshape(self.xglobe.shape), yr.reshape(self.xglobe.shape), zr.reshape(self.xglobe.shape)

        ax = mpl.figure(frameon = False).add_subplot(111,projection = '3d')
        self.globe = ax.plot_surface(xr,yr,zr, facecolors = self.colormap)
        self.globe.remove()
        ax.get_figure().clear()        

class Painter:
    _subjects = {}

    def __init__(self,fig,ax):
        self.board = ax
        self.artists = []     

    def add(cls,objclass,painter):
        Painter._subjects[objclass] = painter

    def get_artists(self):
        return self.artists

    def paint(self, obj, board=None, *args):
        if(board is None):
            board = self.board
        if(isinstance(obj,Orbit)):
            self.paintOrbit(obj,board)
        elif(isinstance(obj,Satellite)):
            self.paintSat(obj,board)
        elif(isinstance(obj,Projectile)):
            self.paintProjectile(obj,board)
        elif(isinstance(obj,Earth)):
            self.paintEarth(obj,board,*args)
        else:
            print(obj)
            print(Orbit)
            for collection in Painter._subjects[obj.__class__](obj,*args).getArtists():
                self.board.add_collection(collection)

    def paintOrbit(self, orbit,board):
        dots = np.linspace(0,2*m.pi,100)
        x_ell = orbit.a*(np.cos(dots)-orbit.e)
        y_ell = orbit.a*m.sqrt(1-orbit.e**2)*np.sin(dots)
        z_ell = np.zeros(x_ell.size)
        orbit_traj = orbit.rel2abs.dot(np.array([x_ell,y_ell,z_ell]))
        board.plot(orbit_traj[0,:],orbit_traj[1,:],orbit_traj[2,:],color = 'g')


    def paintSat(self, sat,board):
        self.satp = satPainter(sat)
        self.artists.extend(self.satp.get_artists())
        board.add_collection(self.satp.body)
        board.add_collection(self.satp.cam)

        # Plot satellite at final position
        r_sat = sat.orbit.getPos(sat.nu)
        v_sat = sat.orbit.getVel(sat.nu)
        m_rot = np.array(sat.attitude.to_matrix())
        # cube_rel = cube()
        # cube_abs = m_rot.dot(cube_rel)
        # xc = sat_rim*cube_abs[0,:].reshape(12,2)+np.ones((12,2))*r_sat[0]
        # yc = sat_rim*cube_abs[1,:].reshape(12,2)+np.ones((12,2))*r_sat[1]
        # zc = sat_rim*cube_abs[2,:].reshape(12,2)+np.ones((12,2))*r_sat[2]
        # poly = board.plot_surface(xc,yc,zc,color = 'm')
        # board.plot_surface(xb,yb,zb,color = 'c')
        # Plot satellite velocity at final position
        mag = 200
        xs,ys,zs = r_sat
        xv,yv,zv = mag*v_sat+r_sat
        varrow = Arrow3D([xs,xv],[ys,yv],[zs,zv], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
        board.add_artist(varrow)

        # Plot satellite Axis
        xd, yd, zd = m_rot.dot(np.array([0,0,sat.width*2]))+r_sat
        sarrow = Arrow3D([xs,xd],[ys,yd],[zs,zd],mutation_scale=20, lw=1, arrowstyle="-|>", color="y")
        board.add_artist(sarrow)

    def paintEarth(self, earth, board, *args):
        globe = False
        if len(args)!=0:
            globe = args[0]
        self.earthp = earthPainter(earth,globe)
        self.artists.extend(self.earthp.get_artists())
        board.add_collection(self.earthp.frame)
        if(globe):
            board.add_collection(self.earthp.globe)
        board.add_artist(self.earthp.meridian)

        # board.plot_wireframe(xg, yg, zg)
        # board.plot(xl,yl,zl,'y')
        # board.plot_surface(xr, yr, zr, rstride=4, cstride=4, facecolors = self.colormap)
        board.set_xlabel('X : vernal point')
        board.set_ylabel('Y : wake dir')
        board.set_zlabel('Z : geo north dir')

    def paintProjectile(self, proj,board):
        board.plot_wireframe(proj.traj[:,0],proj.traj[:,1],proj.traj[:,2],color = 'r')

        # Plot projectile velocity at starting position
        mag = 200
        xs,ys,zs = proj.r_0
        xv,yv,zv = mag*proj.v_0+proj.r_0
        varrow = Arrow3D([xs,xv],[ys,yv],[zs,zv], mutation_scale=20, lw=1, arrowstyle="-|>", color="c")
        board.add_artist(varrow)

    def update(self,index):
        c1 = self.earthp.update(index)
        c2 = self.satp.update(index)
        return c1 + c2
        
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

        

def animate(t,painter,proj,ax, fig, save=False):
    line_p, = ax.plot([], [], [], lw = 2, color = 'r')
    # for collection in painter.get_collections():
    #     ax.add_collection(collection)

    r_proj = proj.traj
    
    def draw(index):
    # update the data
    # t, r_sat, att_sat, r_proj = data
        artists = painter.update(index)    
        line_p.set_data(r_proj[0:index,0],r_proj[0:index,1])
        line_p.set_3d_properties(r_proj[0:index,2])
        return [line_p] + artists

    def init():
        return [line_p] + painter.get_artists()

    ani = animation.FuncAnimation(fig,draw,len(t), blit=False, interval=1,repeat=False, init_func = init)
    if(save):
        ffmpegw = animation.FFMpegWriter()
        f = open('ALE_anim.mp4','w')
        f.close()
        print('Saving animation')
        ani.save(os.path.abspath('ALE_anim.mp4'),ffmpegw)
        print('Animation Saved !')
    mpl.show()
    return ani

def plot(t,earth,sat,proj,**kwargs):
    mpl.close("all")

    animation_on = kwargs.pop('animation',False)
    save_anim = kwargs.pop('save',False)
    see_globe = kwargs.pop('globe',False)

    fig3d = mpl.figure(1)
    ax = fig3d.add_subplot(111, projection='3d')
    ax.set_xlim([-1.1*R_G, 1.1*R_G])
    ax.set_ylim([-1.1*R_G, 1.1*R_G])
    ax.set_zlim([-1.1*R_G, 1.1*R_G])
    renderer = Painter(fig3d,ax)

    renderer.paint(sat.orbit)
    renderer.paint(earth,ax,see_globe)
    renderer.paint(sat)

    
    ani = None
    if(animation_on):
        ani = animate(t,renderer,proj,ax, fig3d, save_anim)
    else:
        renderer.paint(proj,ax)
        
    alt_list = (np.linalg.norm(proj.traj,axis=1)-R_G)*0.001
    
    # Plot descent and sat
    figproj = mpl.figure(2)
    ax = figproj.add_subplot(211)
    ax.plot(t,(np.linalg.norm(proj.traj,axis=1)-R_G)*0.001,'r-')
    ax.plot(t,(np.linalg.norm(sat.traj, axis=1)-R_G)*0.001,'g-')
    mpl.xlabel('Time (s)')
    mpl.ylabel('Altitude (km)')
    mpl.title('Projectile descent in atmosphere')
    
    ax = figproj.add_subplot(212)
    ax.plot(alt_list[1:],earth.atm.rho_hist,'b-')
    mpl.xlabel('Altitude (km)')
    mpl.ylabel('Density (kg/m³)')
    ax.set_yscale('log')
    mpl.title('Density profile during descent')

    # Plot velocity profile
    figvel = mpl.figure(4)
    ax = figvel.add_subplot(211)
    ax.plot((np.linalg.norm(proj.traj,axis=1)-R_G)*0.001,np.linalg.norm(proj.vel,axis=1),'r-')
    mpl.xlabel('Altitude (km)')
    mpl.ylabel('Speed (km/s)')
    mpl.title('Projectile velocity during descent in atmosphere')

    ax = figvel.add_subplot(212)
    ax.plot(proj.Kn,np.linalg.norm(proj.vel,axis=1)[1:],'r-')
    mpl.xlabel('Knudsen (\lambda/d)')
    ax.set_xscale('log')
    mpl.ylabel('Speed (km/s)')
    mpl.title('Projectile velocity at different Knudsen regimes')

    # Plot line of sight
    figlos = mpl.figure(3)
    losvec = figlos.add_subplot(211)
    m_plane = methu.plane_proj(proj.r_0,proj.v_0)
    lostraj = m_plane.dot(sat.lostraj.T).T
    losvec.plot(lostraj[:,0],lostraj[:,1],'c:o')

    losang = figlos.add_subplot(212)
    vec_ref = np.cross(sat.orbit.getVel(sat.nu),-sat.orbit.getPos(sat.nu))
    vec_ref = vec_ref/np.linalg.norm(vec_ref)
    vec_ref = np.tile(vec_ref,sat.traj.shape[0]).reshape(sat.traj.shape)
    ang = methu.angle_vec(-sat.traj,sat.lostraj,vec_ref)
    losang.plot(t,ang,'c')
    mpl.xlabel('Time (s)')
    mpl.ylabel('Angle ( ° )')
    mpl.title('Pointing angle relative to Nadir in orbital dextrogyre reference frame')

    return ani