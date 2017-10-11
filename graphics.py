# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 10:59:16 2017

@author: ASiapan
"""
import numpy as np
from enum import Enum
import os
import PIL
from satclass import *
import matplotlib.pyplot as mpl
import matplotlib
from mpl_toolkits.mplot3d import axes3d, proj3d
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3D
from matplotlib.patches import FancyArrowPatch
from mathutils import Quaternion
import methutil as methu
from methutil import cube
from physenv import *

R_G = 6378 * 10**3 # m 
o = Orbit((30,0,0,0,0))
titleFont = '24'
labelFont = '20'
axesFont = '18'

matplotlib.rc('figure', titlesize = str(int(titleFont)+4))
matplotlib.rc('axes', titlesize = titleFont)
matplotlib.rc('axes', labelsize=labelFont) 
matplotlib.rc('xtick', labelsize=axesFont) 
matplotlib.rc('ytick', labelsize=axesFont)
matplotlib.rc('lines', linewidth = 2.5)
matplotlib.rc('legend', fontsize = labelFont)

class AxName(Enum):
    RHO = "Air Mass Density (kg/m³)"
    TEMPP = "Projectile Temperature (K)"
    TEMPA = "Air Temperature (K)"
    DIAM = "Projectile Diameter (cm) "
    Ma = "Mach Number"
    Kn = "Knudsen Number"
    Kn_mix = "Mixed Knudsen Number"
    Re = "Reynolds Number"
    tr = "(Projectile)/(Air)\n Temperature Ratio"
    ALT = "Altitude (km)"
    SPEED = "Speed (km/s)"
    C_D = "Drag Coefficient"
    NINC = "Incident gas flow\n on surface (#/s) "
    sigma_d = "Surface Accomodation Coefficient"
    HEAT = "Total Integrated Heat Flow (W)"
    HEATC = "Collisional Heat Flow (W)"
    HEATR = "Radiative Heat Flow (W)"

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
        rotmat = np.array(earth.tilt_rot.to_matrix())
        xgr, ygr, zgr = rotmat.dot(np.vstack((xg.flatten(),yg.flatten(),zg.flatten())))
        xg = xgr.reshape(xg.shape)
        yg = ygr.reshape(yg.shape)
        zg = zgr.reshape(zg.shape)



        ax = mpl.figure(frameon = False).add_subplot(111,projection = '3d')
        self.frame = ax.plot_wireframe(xg,yg,zg)
        self.frame.remove()
        ax.get_figure().clear()
        

        # self.phi0 = earth.phi[0]
        self.rotation = Quaternion(earth.axis,earth.rotang)
        # self.rotation.angle = self.phi0
        m_rot = np.array((earth.tilt_rot*Quaternion([0,0,1],earth.rotang)).to_matrix())
        # self.xm = R_G * np.cos(self.phi0)*np.sin(theta)
        # self.ym = R_G * np.sin(self.phi0)*np.sin(theta)
        self.xm = R_G * np.sin(theta)
        self.zm = R_G * np.cos(theta)

        self.xm, self.ym, self.zm = m_rot.dot(np.array((self.xm,np.zeros((self.xm.shape)),self.zm)))

        self.meridian = Line3D(self.xm,self.ym,self.zm,color = 'y')
        self._see_globe = globe
        if(globe):
            self.load_globe()


    def update(self,index):
        self.rotation.angle = self.earth.phi[index]-self.earth.rotang
        m_rot = np.array((self.rotation).to_matrix())
        
        if(self._see_globe):
            m_rot = np.array(self.rotation.to_matrix())
            tmp = np.vstack((self.xglobe.reshape((1,self.xglobe.size)), self.yglobe.reshape((1,self.yglobe.size)), self.zglobe.reshape((1,self.zglobe.size))))
            xr,yr,zr = m_rot.dot(tmp)
            xr,yr,zr = xr.reshape(self.xglobe.shape), yr.reshape(self.xglobe.shape), zr.reshape(self.xglobe.shape)
            self.globe.set_verts(convert2vert(xr,yr,zr))

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
        self.globe = ax.plot_surface(xr,yr,zr, rstride=4, cstride=4, facecolors = self.colormap)
        self.globe.remove()
        self.globe.set_zorder(0)
        # ax.get_figure().clear()
        mpl.close(ax.get_figure())
        print("hello")   

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
        if(isinstance(obj,Orbit) or obj.__class__.__name__=='Orbit'):
            self.paintOrbit(obj,board,*args)
        elif(isinstance(obj,Satellite) or obj.__class__.__name__=='Satellite'):
            self.paintSat(obj,board,*args)
        elif(isinstance(obj,Projectile) or obj.__class__.__name__=='Projectile' or obj.__class__.__name__=='Ghostectile'):
            self.paintProjectile(obj,board,*args)
        elif(isinstance(obj,Earth) or obj.__class__.__name__=='Earth'):
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
        board.plot_wireframe(orbit_traj[0,:],orbit_traj[1,:],orbit_traj[2,:],color = 'g')

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
        if(globe):
            board.add_collection(self.earthp.globe)
            print('globe ',self.earthp.globe, " ",self.earthp.globe.zorder)
        else:
            board.add_collection(self.earthp.frame)
        board.add_artist(self.earthp.meridian)

        board.set_xlabel('X : vernal point')
        board.set_ylabel('Y : wake dir')
        board.set_zlabel('Z : geo north dir')

    def paintProjectile(self, proj,board, c = 'r'):
        board.plot_wireframe(proj.traj[:,0],proj.traj[:,1],proj.traj[:,2],color = c)

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

def palette(color_0,n_colors, natural = False):
    if(natural): 
        colors = []
        for i in range(n_colors):
            colors.append('C'+str(i))
    else:
        increment  = 1/n_colors
        mask = np.array((-0.2,1,0.25))
        colors = [color_0]
        for c in range(1,n_colors):
            color = color_0+c*increment*mask
            color = color%1
            colors.append(color)
    return colors

def markplot(ax,x,y,id = None, marker = 'x',color = 'r'):
    if(not id): id = range(len(x))
    xi = [x[i] for i in id]
    yi = [y[i] for i in id]
    artist = ax.scatter(xi,yi, marker = marker, s = 100, c = color)
    return artist

def unitplot(ax,AxX,AxY,title,legend=None, padding = 0,color = 'r'):
    ax.plot(AxX[0],AxY[0],color = color)
    ax.set_xlabel(AxX[1],labelpad = padding)
    ax.set_ylabel(AxY[1])
    ax.set_title(title)
    

def plot(t,earth,sat,proj,hypbox,**kwargs):
    mpl.close("all")
    figind = 1
    axes = {}


    animation_on = kwargs.pop('animation',False)
    save_anim = kwargs.pop('save',False)
    see_globe = kwargs.pop('globe',False)
    imark, marks = kwargs.pop('marks',(0,[]))
    db = kwargs.pop('dsmcdata',[])


    fig3d = mpl.figure(figind)
    ax = fig3d.add_subplot(111, projection='3d')
    ax.set_xlim([-1.1*R_G, 1.1*R_G])
    ax.set_ylim([-1.1*R_G, 1.1*R_G])
    ax.set_zlim([-1.1*R_G, 1.1*R_G])
    renderer = Painter(fig3d,ax)
    axes["Animation"] = ax

    renderer.paint(sat.orbit)
    renderer.paint(earth,ax,see_globe)
    renderer.paint(sat)

    n_traj = len(proj)
    colors = kwargs.pop('palette',palette(np.array((1.,0,0)),n_traj))
    
    ani = None
    if(animation_on):
        ani = animate(t,renderer,proj[0],ax, fig3d, save_anim)
    else:
        for p in range(n_traj):
            renderer.paint(proj[p],ax,colors[p])
    figind += 1

    name_list = []
    alt_mark = [r.alt for r in db]
    C_Dmark = [r.C_D for r in db]
    qmark = [r.q for r in db]

    for p in range(n_traj):
        name_list.append(proj[p].id)
    
    # Plot descent and sat
    figproj = mpl.figure(figind)
    figind += 1
    ax = figproj.add_subplot(211)
    for p in range(n_traj): ax.plot(proj[p].time,hypbox[p].altline,color = colors[p])
    if marks : markplot(ax,proj[imark].time,hypbox[imark].altline,marks,color=colors[imark])
    ax.plot(t,(np.linalg.norm(sat.traj, axis=1)-R_G)*0.001,'g-')
    ax.set_xlabel('Time (s)',labelpad = -50)
    ax.set_ylabel('Altitude (km)')
    ax.set_title('Projectile descent in atmosphere')
    axes["Descent"] = ax
    
    ax = figproj.add_subplot(212)
    for p in range(n_traj): ax.plot(hypbox[p].altline[1:],hypbox[p].rholine,color = colors[p])
    if marks : markplot(ax,hypbox[imark].altline[1:],hypbox[imark].rholine,marks,color=colors[imark])
    ax.set_xlabel('Altitude (km)')
    ax.set_ylabel('Density (kg/m³)')
    ax.set_yscale('log')
    ax.set_title('Density profile during descent')
    axes["Density"] = ax

    # Plot velocity profile
    figvel = mpl.figure(figind)
    figind += 1
    ax = figvel.add_subplot(211)
    for p in range(n_traj): ax.plot(hypbox[p].altline,hypbox[p].speedline, color = colors[p])
    if marks : markplot(ax,hypbox[imark].altline,hypbox[imark].speedline,marks,color = colors[imark])
    ax.set_xlabel('Altitude (km)',labelpad = -50)
    ax.set_ylabel('Speed (km/s)')
    ax.set_title('Projectile velocity during descent in atmosphere')
    if(n_traj==1):
        ax2 = ax.twinx()
        for p in range(n_traj): ax2.plot(hypbox[p].altline[1:],hypbox[p].machline, color = 'b',linestyle = "dashed")
        for p in range(n_traj): ax.plot([],[], color = 'b',linestyle = "dashed")
        ax2.set_ylabel(AxName.Ma.value)
        ax.legend(['Velocity', 'Mach number'])
    else: ax.legend([str(p) for p in proj])
    axes["FlightPath"] = ax

    ax = figvel.add_subplot(212)
    for p in range(n_traj): ax.plot(hypbox[p].Knline,hypbox[p].speedline[1:], color = colors[p])
    ax.set_xlabel(AxName.Kn.value)
    ax.set_xscale('log')
    ax.set_ylabel('Speed (km/s)')
    ax.set_title('Projectile velocity at different Knudsen regimes')
    if(n_traj>1): ax.legend([str(p) for p in proj])
    axes["VelKn"] = ax

    # Plot line of sight
    figlos = mpl.figure(figind)
    figind += 1
    losvec = figlos.add_subplot(211)
    # for p in range(n_traj):
    #     if(proj[p].__class__.__name__ != 'Ghostectile'):
    #         m_plane = methu.plane_proj(proj[p].r_0,proj[p].v_0)
    #         lostraj = m_plane.dot(sat.lostraj[proj[p]].T).T
    #         losvec.plot(lostraj[:,0],lostraj[:,1], color = colors[p])
    # losvec.set_xlabel('Relative X')
    # losvec.set_ylabel('Relative Y')
    # losvec.set_title('Pointing vector direction in the projectile initial plane reference frame')

    losang = figlos.add_subplot(212)
    vec_ref = np.cross(sat.orbit.getVel(sat.nu),-sat.orbit.getPos(sat.nu))
    vec_ref = vec_ref/np.linalg.norm(vec_ref)
    vec_ref = np.tile(vec_ref,sat.traj.shape[0]).reshape(sat.traj.shape)
    for p in range(n_traj):
        if(proj[p].__class__.__name__ != 'Ghostectile'):
            ang = methu.angle_vec(-sat.traj,sat.lostraj[proj[p]],vec_ref)*180/np.pi
            losang.plot(t,ang,color = colors[p])
    losang.set_xlabel('Time (s)')
    losang.set_ylabel('Angle ( ° )')
    losang.set_title('Pointing angle relative to Nadir in orbital dextrogyre reference frame')
    if(n_traj>1): losang.legend([str(p) for p in proj])
    axes["PointingAngle"] = losang

    # Plot drag coefficient
    figdrag = mpl.figure(figind)
    figind += 1
    ax = figdrag.add_subplot(211)
    for p in range(n_traj): ax.plot(hypbox[p].altline[1:],hypbox[p].dragline,color = colors[p])
    markplot(ax,alt_mark,C_Dmark,marker='o',color='g')
    ax.set_xlabel(AxName.ALT.value)
    ax.set_ylabel(AxName.C_D.value)
    ax.set_title('Projectile drag coefficient at different altitudes')
    if(n_traj>1): ax.legend([str(p) for p in proj]+["DSMC instances for Ghost"])
    axes["C_D"] = ax

    # Plot temperature and heating
    figheat = mpl.figure(figind)
    figind += 1
    ax = figheat.add_subplot(211)
    for p in range(n_traj): ax.plot(hypbox[p].altline[:-2],proj[p].temp[:-2],color = colors[p])
    ax.set_xlabel('Altitude (km)')
    ax.set_ylabel('Temperature (K) _ solid _')
    #ax2 = figheat.add_subplot(212, sharex = ax)
    ax2 = ax.twinx()
    for p in range(n_traj): ax2.plot(hypbox[p].altline[1:], hypbox[p].heatline, color = colors[p], linestyle = '--')
    ax2.set_ylabel('Heat Flow (W) -- dashed --')
    markplot(ax2,alt_mark,qmark,range(len(alt_mark)),'o','g')
    ax.set_title('Projectile temperature and heat at different altitudes')
    if(n_traj>1): ax2.legend([str(p) for p in proj])#+["DSMC instances for Ghost"])
    axes["T"] = ax
    axes["q"] = ax2

    # Plot incoming collisions and ablation
    figcoll = mpl.figure(figind)
    figind += 1
    ax = figcoll.add_subplot(211)
    for p in range(n_traj): ax.plot(hypbox[p].altline[:],np.array(proj[p].diam)*100,color = colors[p])
    ax.set_xlabel(AxName.ALT.value, labelpad = -50)
    ax.set_ylabel(AxName.DIAM.value)
    ax.set_title('Projectile diameter along altitude')
    if(n_traj>1): ax.legend([str(p) for p in proj])
    axes["Diam"] = ax
    ax = figcoll.add_subplot(212, sharex = ax)
    for p in range(n_traj): ax.plot(hypbox[p].altline[1:], hypbox[p].collline, color = colors[p])
    ax.set_xlabel(AxName.ALT.value)
    ax.set_ylabel(AxName.NINC.value)
    ax.set_title('Oncoming particle flux on the projectile\'s surface')
    if(n_traj>1): ax.legend([str(p) for p in proj])
    axes["Ni"] = ax

    # Plot regimes - Kn, Re, Ma and t_r
    figreg = mpl.figure(figind)
    figind += 1
    figreg.suptitle('Similarity Parameters along altitude')
    ax = figreg.add_subplot(221)
    els = []
    for p in range(n_traj): els.extend(ax.plot(hypbox[p].altline[1:],hypbox[p].machline,color = colors[p]))
    if marks : els.append(markplot(ax,hypbox[imark].altline,hypbox[imark].machline,marks,color=colors[imark]))
    ax.set_xlabel(AxName.ALT.value)
    ax.set_ylabel(AxName.Ma.value)
    # if(n_traj>1): ax.legend([str(p) for p in proj])
    axes["Ma"] = ax

    ax = figreg.add_subplot(222, sharex = ax)
    for p in range(n_traj): ax.plot(hypbox[p].altline[1:],hypbox[p].Reline,color = colors[p])
    if marks : markplot(ax,hypbox[imark].altline,hypbox[imark].Reline,marks,color=colors[imark])
    ax.set_xlabel(AxName.ALT.value)
    ax.set_ylabel(AxName.Re.value)
    # if(n_traj>1): ax.legend([str(p) for p in proj])
    axes["Re"] = ax

    ax = figreg.add_subplot(223, sharex = ax)
    for p in range(n_traj): ax.plot(hypbox[p].altline[1:],hypbox[p].Knline,color = colors[p])
    if marks : markplot(ax,hypbox[imark].altline,hypbox[imark].Knline,marks,color=colors[imark])
    ax.set_xlabel(AxName.ALT.value)
    ax.set_ylabel(AxName.Kn.value)
    ax.set_yscale('log')
    # if(n_traj>1): ax.legend([str(p) for p in proj])
    axes["Kn"] = ax

    ax = figreg.add_subplot(224, sharex = ax)
    for p in range(n_traj): ax.plot(hypbox[p].altline[1:],hypbox[p].trline,color = colors[p])
    if marks : markplot(ax,hypbox[imark].altline,hypbox[imark].trline,marks,color=colors[imark])
    ax.set_xlabel(AxName.ALT.value)
    ax.set_ylabel(AxName.tr.value)
    # if(n_traj>1): ax.legend([str(p) for p in proj])
    axes["tr"] = ax
    if(n_traj>1): figreg.legend(els,[str(p) for p in proj])
    else: figreg.legend(els,['Projectile','Ghostectile','DSMC instances'])

    return ani, axes

def quickplot(*args,same = False, colors = ['r']):
    assert len(args)%3==0
    nplots = int(len(args)/3)
    if (len(colors)!= nplots) : colors *= nplots
    if(same==True): 
        fig, axs = mpl.subplots(2,1)
        axs = tuple([axs[0]]*nplots)
    else : fig, axs = mpl.subplots(nplots,1)
    if(nplots == 1) : axs = (axs,)
    for pl in range(nplots):
        AxX = args[3*pl+0]
        AxY = args[3*pl+1]
        title = args[3*pl+2]
        if(all(type(y)==list for y in AxY[0])):
            for y in AxY[0]:
                axs[pl].plot(AxX[0],y,color = colors[pl])
        else: axs[pl].plot(AxX[0],AxY[0],color = colors[pl])
        axs[pl].set_xlabel(AxX[1])#,labelpad = -50)
        axs[pl].set_ylabel(AxY[1])
        if(len(AxY)>2): axs[pl].set_yscale(AxY[2])
        if(len(AxX)>2): axs[pl].set_xscale(AxX[2])
        axs[pl].set_title(title)
    return axs

def compare(sat,proj,hypbox):
    n_traj = len(proj)
    t = hypbox[0].timeline
    colors = palette(np.array((1.,0,0)),n_traj)
    ang_list = []
    for p in range(n_traj):
        if(proj[p].__class__.__name__ != 'Ghostectile'):
            vec_ref = np.cross(sat.orbit.getVel(sat.nu),-sat.orbit.getPos(sat.nu))
            vec_ref = vec_ref/np.linalg.norm(vec_ref)
            vec_ref = np.tile(vec_ref,sat.lostraj[proj[p]].shape[0]).reshape(sat.lostraj[proj[p]].shape)
            ang_list.append(methu.angle_vec(-sat.traj,sat.lostraj[proj[p]],vec_ref)*180/np.pi)
    ang_comp = ang_list[0]-ang_list[1]
    alt_comp = (hypbox[0].altline-hypbox[1].altline)/hypbox[0].altline*100
    pos_comp = np.linalg.norm(proj[0].traj-proj[1].traj, axis = 1)/np.linalg.norm(proj[0].traj,axis = 1)*100
    vel_comp = np.linalg.norm(proj[0].vel-proj[1].vel, axis = 1)/hypbox[0].speedline*100
    
    # Plot drag coefficient
    # figdrag = mpl.figure(5)
    fig, ax = mpl.subplots()#figdrag.add_subplot(311)
    for p in range(n_traj): ax.plot(hypbox[p].altline[1:],hypbox[p].dragline,color = colors[p])
    ax.set_xlabel('Altitude (km)')
    ax.set_ylabel('Drag coefficient')
    ax.legend(['C_Drag = 1','C_Drag = 2'])
    ax.set_title('Projectile drag coefficient at different altitudes')
    
    fig, ax = mpl.subplots(2,1)#figdrag.add_subplot(312)
    fig.suptitle('Comparison between two trajectories with different drag coefficients (C_D = 1 & 2)',fontsize = 20)
    ax[0].plot(t,alt_comp,'r-')
    ax[0].set_xlabel('Time (s)', fontsize = labelFont)
    ax[0].set_ylabel('Altitude relative \n difference (%)', fontsize = labelFont)
    ax[1].plot(hypbox[0].altline*0.5+hypbox[1].altline*0.5,pos_comp,'r-')
    ax[1].set_xlabel('Altitude (mean) (km)', fontsize = labelFont)
    ax[1].set_ylabel('Position relative \n difference (%)', fontsize = labelFont)
    
    fig, ax = mpl.subplots()#figdrag.add_subplot(313)
    ax.plot(t,vel_comp,'r-')
    ax.set_xlabel('Time (s)', fontsize = labelFont)
    ax.set_ylabel('Velocity relative \n difference (%)', fontsize = labelFont)
    ax.get_xticklabels

    fig, ax = mpl.subplots()
    ax.plot(t,ang_comp,'r-')
    ax.set_xlabel('Time (s)', fontsize = labelFont)
    ax.set_ylabel('Angle absolute \n difference ( ° )', fontsize = labelFont)
    ax.set_title('Pointing angle relative to Nadir in orbital dextrogyre reference frame', fontsize = titleFont)    

def plot_atm(atmosphere,h_range = np.linspace(50,400,801)):
    lat, lon, name = Reference.get_refloc()
    h_scale, rho_scale, T_scale = atmosphere.profile(h_range,lat, lon)
    mpl.figure()
    mpl.plot(h_scale,rho_scale,'b-')
    ax.set_title('Density profile at {0:+.2f} Latitude, {1:+.2f} Longitude ({2})'.format(lat*180/m.pi,lon*180/m.pi,name))
    ax.set_xlabel('Altitude (km)')
    ax.set_ylabel('Mass Density (kg/m³)')