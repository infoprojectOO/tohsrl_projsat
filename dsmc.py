# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:26:28 2017

@author: ASiapan
"""

# DSMC code implementation

import math as m
import numpy as np
import random as rand
import methutil as methu
from itertools import product
from dsmcgraph import showfield
import traceback
import os, sys
import os, pathlib, shlex
import subprocess
import re
import time
import param as cst
from enum import Enum
import pathlib
import glob
import sys,os
from queue import Queue
from datetime import datetime
from threading import Thread
import pickle


path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from olog import olog
from dump import dump


sparta_folder = "\\dsmc\\"

class SimInput(Enum):
    '''Keys to replace DSMC input parameters with the desired values'''
    ID = "id", "0"
    AIR_NDENSITY = "nv_air", "1e19"
    AIR_TEMPERATURE = "T_air", "200"
    SPEED = "v_air", "7000"
    PROJ_TEMPERATURE = "T_proj" 
    PROJ_SURFSTATE = "sigma_d", "1"
    PROJ_SIZE = "d_proj", str(cst.d_p)
    TIMESAMP = "t_ref", r"${t_ref}", r"${t_step}"
    PROJ_CALCAP = "Ccal_proj"

    def __init__(self,value,default=None, fine=None):
        self._value = value
        self.default = default
        self.fine = fine
    
    @property
    def value(self):
        return self._value

class SimOutput(Enum):
    '''Keys to specify the output format for the DSMC runs'''
    SURF_FILE = "surfout" # Surface vectors (p,T,heat, etc.)
    FLOW_FILE = "flowout" # Flow grid arrays (p,T,rho, etc.)
    LOG_FILE = "logfile" # Macro values (Total drag, Total heat, Total impact flux, etc.)
    SURF_OUT = "surf_comp" # what to output for the surface
    FLOW_OUT = "flow_comp" # what to output for the flow

class Index:
    def __init__(self,value):
        self.value = value

class SimOption(Enum):
    '''Options that specify the outputs to be generated'''
    DRAG_COEFF = ("fx","fx",SimOutput.LOG_FILE,SimOutput.SURF_OUT)
    HEAT = ("qtot","qtot",SimOutput.LOG_FILE,SimOutput.SURF_OUT)
    SCOLL = ("Ni","nwt",SimOutput.LOG_FILE,SimOutput.SURF_OUT)
    AIR_DENSITY = ("n_flow","n", SimOutput.FLOW_FILE,SimOutput.FLOW_OUT)
    AIR_TEMPERATURE = ("T_flow","temp",SimOutput.FLOW_FILE,SimOutput.FLOW_OUT)
    AIR_TEMPROT = ("Trot_flow","trot",SimOutput.FLOW_FILE,SimOutput.FLOW_OUT)
    AIR_TEMPVIB =  ("Tvib_flow","tvib",SimOutput.FLOW_FILE,SimOutput.FLOW_OUT)   
    
    def __init__(self,value,subst,*deps):
        self._file = deps[0]
        self._value = value
        self.key = subst
        self.index = Index("i"+value)
        self.dependencies = deps
        
    @property
    def value(self):
        return self._value
    @property
    def fileid(self):
        return self._file

    def __add__(self,so):
        return self.value + so.value

class SimID:
    '''Generate automatic DSMC simulation naming format for files used'''

    def __init__(self,name):
        self.key = name

    def __add__(self,other):
        return SimID(self.key+','+other.key)

    def __str__(self):
        return self.key

    class keys(Enum):
        alt = "{alt:.2f}km"
        Ma = "Ma_{Ma:d}"
        Kn = "Kn_{Kn:d}"
        Re = "Re_{Re:d}"
        num = "{0:d}-"
        proj = "{proj:s}"
        ID = "id-{id:s}"
        tw = "tw_{tw:.2f}"

        def __init__(self,name):
            self.key = name

        def __str__(self):
            return self.key

        def __add__(self,other):
            if(self.key[-1]!="-"):
                return SimID(self.key+other.key)
            else:
                return SimID(self.key+','+other.key)

class DSMC:
    '''Manages DSMC simulations, stores them in a database'''
    dbfile = "datab"
    dbformat = str(SimID.keys.num+SimID.keys.Ma+SimID.keys.Kn+SimID.keys.Re+SimID.keys.tw)
    simformat = str(SimID.keys.num+SimID.keys.alt+SimID.keys.Ma+SimID.keys.ID)
    stamp = "ver:0"#methu.version_gen("ver:",2)
    dsmcID = 0

    def __init__(self, datafile = dbfile, idformat=dbformat, id = None):
        if(id): self.ID = id
        else: self.ID = str(DSMC.dsmcID)
        DSMC.dsmcID += 1
        self.path2dsmc = os.getcwd()+sparta_folder
        self.runner = Simulator(self.pull)
        self.stop = False
        # self.sim_options = [SimOption.DRAG_COEFF,SimOption.HEAT,SimOption.SCOLL]
        self.regimes = {}
        self.id2key = {}
        self.simlist = []
        self.datafile = datafile+self.ID
        self.format = idformat
        # self.running = {}
        # self.finished = {}
        # self.simind = 0
        # self.io = SpartaIO()
        # self.io.define_simulation(self.sim_options)
        # self.sentinel = Thread(target = self.watch)
        self.loadData(self.datafile)
        # self.sentinel.start()
        self.nsim = 0

    def loadData(self,fileid):
        '''Reads the regime data from the database and loads it into the virtual memory

        Arguments
        ---------
        fileid : string - full name of the database file

        Returns
        -------
        AeroRegime{...} - regime dictionary with key according to simulation name formatting
        '''
        datafile = pathlib.Path('.'+sparta_folder+fileid)
        if(datafile.exists()):
            try:
                with open(str(datafile),'rb') as f:
                    self.database = pickle.load(f)
            except:
                print('Couldn\'t open previous database, creating a new one')
                self.database = {}
        else:
            self.database = {}
            datafile.touch(exist_ok = False) # Do not allow to overwrite an existing database, abort immediately
        for key in self.database.copy().keys():
            if(self.database[key].ver != self.stamp):
                self.database.pop(key)    

    def add_sim(self,aeroreg,t_step,force=False):
        '''Add a simulation instances onto the stack of orders

        Arguments
        ---------
        aeroreg : AeroRegime object - represents the regime instance to be simulated (input parameters)
        t_step : specifies the sampling time step, relates to the trade-off between accuracy and rapidity
        --------- Keyword args
        force : boolean - Should we overwrite the simulation if a similar one is already found present in the database ?
        '''
        self.nsim += 1
        aerokey = self.format.format(self.nsim,alt = aeroreg.alt,Ma=int(np.round(aeroreg.Ma)),Kn=int(np.round(aeroreg.Kn)),Re = int(np.round(aeroreg.Re)), tw = aeroreg.tr, id= self.ID)
        if(not force and aerokey in self.database and self.database[aerokey].ver == self.stamp): return False
        aeroreg.ver = "None" # The regime has not yet been numerically simulated, it does not contain any simulation results yet
        fresh = len(self.simlist)==0 #Is this the first simulation to the processed ?
        self.simlist.append(aerokey)
        self.database[aerokey] = aeroreg
        simid = self.runner.order_sim(aeroreg,fresh=fresh)
#        elif(len(self.finished.values()) == 1):
#            self.runner.order_sim(aeroreg,simid=aerokey)
            # self.launch(simid,[self.sparta_exe,"-in",simname])
        self.id2key[simid] = aerokey
        return True

    def pull(self,simid,aeroreg):
        ''' Retrieve data from simulation results and stamp the regime that was processed with a version number

        Arguments
        ---------
        simid : string - how the simulation was named
        aeroreg : AeroRegime object - the regime that carries the simulation results information
        '''
        aeroreg.ver = self.stamp
        aerokey = self.id2key[simid]
        self.database[aerokey] = aeroreg

    def exertion(self,aeroreg):
        ''' Computes the drag, heat and impact flux requested by the projectile trajectory prediction algorithm'''
        alt = aeroreg.alt
        keyids = self.find_closest(aeroreg)
        ready = False
        first = True
        if all([id in self.database for id in keyids]):
            ready = all([self.database[id].ver == self.stamp for id in keyids])
        while (not (ready) or len(self.database)==0) and not self.stop:
            if(first): 
                first = False
                print("Waiting for simulation at {:d} km".format(int(alt)),end = '') 
            else: print('.',end = '')
            time.sleep(3) # waiting for simulation
            keyids = self.find_closest(aeroreg)
            if all([id in self.database for id in keyids]):
                ready = all([self.database[id].ver == self.stamp for id in keyids])
        return self.combine(aeroreg,[self.database[id] for id in keyids])

    def combine(self,reg,basis):
        ''' Interpolates the current regime values with a set of DSMC simulated regime instances

        Arguments
        ---------
        reg : AeroRegime object - the current regime
        basis : AeroRegime[] list - the regimes served as nodes for interpolation

        Returns
        -------
        drag coefficient, heat coefficient and impact flux coefficient
        '''
        weight = np.array([])
        drags = []
        heats = []
        incs = []
        for r in basis:
            dist = reg.distance(r)
            weight = np.hstack((weight,1/dist))
            drags.append(r.C_D)
            heats.append(r.C_H)
            incs.append(r.C_Ni)
        tot = sum(weight)
        weight = weight/tot
        C_D = np.array(drags).dot(weight)
        heat_C = np.array(heats).dot(weight)
        onflow_C = np.array(incs).dot(weight)
        return C_D, heat_C, onflow_C

    def find_closest(self,reg,n=4):
        ''' Find closest regimes in database to the regime given in parameter

        Arguments
        ---------
        reg : AeroRegime object : the regime that seeks for neighbours
        --------- Keyword args
        n : number of neighbouring regimes sought
        '''
        dist = {}
        for id,r in self.database.items():
            d = reg.distance(r)
            dist[id] = d
        try:
            ids = methu.minima(dist,n)
        except Exception as e:
            #print(e)
            ids = dist.keys() # The entire database is used if a problem occurs (the only solution left)
        return ids

    def abort_all(self):
        self.stop = True
        self.runner.abort_all()

    def save_data(self):
        ''' Overwrites the database file with its updated database'''
        with open(self.path2dsmc+self.datafile,'wb') as f:
            pickle.dump(self.database,f)

class Simulator:
    ''' Handles all DSMC simulations, file naming, simulation queuing, retrieving data outputs'''

    MAX_SIM_PARA = 1
    sim_folder = "dsmc"
    exefile = os.getcwd()+sparta_folder+"spa_radial"
    batch_file_name = "dsmc_cmd.bat"
    cygwin_path = "C:\\cygwin\\bin\\mintty.exe"
    bash_path = "C:\\cygwin\\bin\\bash.exe"
    simformat = str(SimID.keys.num+SimID.keys.alt+SimID.keys.Ma+SimID.keys.ID)
    dsmcID = 0 # Used for indexing simulations

    def __init__(self, callback = None):
        self.sim_options = [SimOption.DRAG_COEFF,SimOption.HEAT,SimOption.SCOLL]#,SimOption.AIR_DENSITY,SimOption.AIR_TEMPERATURE]
        self.io = SpartaIO()
        self.io.define_simulation(self.sim_options)
        self.running = {}
        self.queue = Queue()
        self.simlist = []
        self.finished = {}
        self.outputs = {}
        self.regimes = {}
        self.sentinel = Thread(target = self.watch)
        self.isim = 1
        self.callback = callback
        self.ID = str(Simulator.dsmcID)
        Simulator.dsmcID += 1

    def reset(self):
        '''Clear simulation list, results, abort all current simulations and restart all threads'''
        self.abort_all()
        self.sentinel = Thread(target = self.watch)
        self.simlist.clear()
        self.running.clear()
        self.finished.clear()
        self.queue = Queue()
        self.isim = 1

    def abort_all(self):
        ''' Terminate immediately all currently running simulations '''
        self.watching = False
        for f,p in self.running.items():
            p.terminate()       

    def replay(self,nsim):
        ''' Replay any simulation according to its numbering index

        Arguments
        ---------
        nsim : number index of simulation in the projectile descent course

        Returns
        -------
        Subprocess object - the simulation rerunning process
        '''
        simfile = self.simlist[nsim]["simfile"]
        simname = self.simlist[nsim]["simname"]
        instr = [self.sparta_exe,"-in",simfile]
        sim = subprocess.Popen(instr,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = self.sparta_folder, universal_newlines = True, shell=False)
        return sim

    def set_trace(self,recbox):
        ''' Used for recording all regimes used in the simulation inputs'''
        self.recbox = recbox
        self.thermovody = recbox.movody;

    def order_sim(self,aeroreg,t_step = SimInput.TIMESAMP.default, simid = None, fresh=False):
        ''' Add a simulation order to the queue of simulations

        Arguments
        ---------
        aeroreg : AeroRegime object - regime conditions for input conditions
        ---------Keyword args
        t_step : sampling time step < default defined in script
        simid : string - unique simulation identification string < default uses standard format naming
        fresh : boolean - Should the simulation be launched serially for its successors ? < default is parallel
        '''
        if simid is None : simid = self.simformat.format(self.isim,alt = aeroreg.alt,Ma=int(np.round(aeroreg.Ma)),Kn=int(np.round(aeroreg.Kn)),Re = int(np.round(aeroreg.Re)), tw = aeroreg.tr, id= self.ID)        
        sim_params = self.buildparams(simid,aeroreg,t_step)
        simname = self.io.write(sim_params,fresh)
        self.queue.put(sim_params)
        self.simlist.append(sim_params)
        print("Preparing Simulation : ",simname.split("\\")[-1])
        # simid = sim_params["simname"]
        self.regimes[simid] = aeroreg
        if(len(self.running.values()) < self.MAX_SIM_PARA):
            simparams = self.queue.get()
            self.launch(simid,[self.exefile,"-in",simname])
        self.isim += 1
        return simid

    def buildparams(self,key,aeroreg,t_step = SimInput.TIMESAMP.default):
        '''Extract parameters from regime and translate in input parameter bundle

        Arguments
        ---------
        key : SimInput enum - unique parameter input key
        aeroreg : AeroRegime object - contains values for each key
        --------- Keyword args
        t_stap : DSMC sampling time step
        '''
        alt = aeroreg.alt
        Ma = aeroreg.Ma
        param = {SimInput.ID:key}
        param[SimInput.AIR_NDENSITY] = aeroreg.nv
        param[SimInput.AIR_TEMPERATURE] = aeroreg.T
        param[SimInput.SPEED] = aeroreg.U
        param[SimInput.PROJ_TEMPERATURE] = aeroreg.Tr
        param[SimInput.PROJ_SIZE] = aeroreg.L
        param[SimInput.PROJ_SURFSTATE] = aeroreg.sigma_d
        param[SimInput.PROJ_CALCAP] = aeroreg.C_cal*aeroreg.m
        param[SimInput.TIMESAMP] = t_step
        return param

    def watch(self):
        ''' Start watching over DSMC simulation threads, waiting to extract any output after they finish their run'''
        self.watching = True
        self.delay = 5 # s - used for sleeping time 
        while self.watching:
            for fileid,proc in self.running.copy().items(): # Look over running processes
                if proc.poll(): # Process finished with or without error
                    print(fileid," Simulation results : ",proc.communicate()[1])
                    self.running.pop(fileid)
                    self.retrieve_results(fileid)
                    self.next_run()
                    self.finished[fileid] = proc # Keep record of process threads launched
                    if(len(self.running)==0 and self.queue.qsize()==0) : self.watching = False # No one else is running and no one else is waiting to be launched
                else: # Process is not finished
                    try:
                        proc.communicate(timeout = self.delay) # Try to extract any information
                    except subprocess.TimeoutExpired: # Process does not respond
                        continue
            # time.sleep(self.delay)

    def retrieve_results(self,fileid):
        ''' Extract simulation results and write them into the original input regime object

        Arguments
        ---------
        fileid : string - simulation unique identifier
        '''
        try:
            print("Retrieving from ",fileid)
            output = self.io.read(fileid,self.sim_options)
            reg = self.regimes[fileid]
            reg.C_D = np.mean(output[SimOption.DRAG_COEFF][1:]) # Mean value over whole simulation samples (excluding the first row which is initialisation)
            reg.C_H = np.mean(output[SimOption.HEAT][1:])/(0.5*reg.air.rho*reg.S*reg.U**3)#(reg.air.p*reg.S*reg.air.v_therm)
            reg.C_Ni = np.mean(output[SimOption.SCOLL][1:])/(0.5*reg.nv*reg.S*reg.U)#/(reg.nv*reg.S*reg.air.v_therm)
            self.outputs[fileid] = output
            if self.callback is not None: self.callback(fileid,reg)
        except:
            traceback.print_exc()
            self.abort_all()

    def next_run(self):
        ''' Launch next simulation(s) waiting in the queue depending on the amount of parallel batches available'''
        while(len(self.running) < self.MAX_SIM_PARA and self.queue.qsize()>0):
            simparam = self.queue.get()
            simfile = simparam["simfile"] # Input file for the simulation
            simname = simparam["simname"] # Simulation identifier
            self.launch(simname,[self.exefile,"-in",simfile])
        # print("launched next : returning to watch   ", self.running, "  queue size :   ",self.queue.qsize())

    def launch(self,simid,instr):
        ''' Launch a process

        Arguments
        ---------
        simid : string - simulation identifier
        instr : string - instruction given to terminal
        '''
        if(not self.sentinel.is_alive()): self.sentinel.start()
        print("Launching Simulation : "+simid)
        sim = subprocess.Popen(instr,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = self.sim_folder, universal_newlines = True, shell=False)
        self.running[simid] = sim
        return sim

    def prepare_batch_file(self,filename):
        f = pathlib.Path(filename)
        instr = "spa_serial < in.axinit"
        # if(not f.exists()):
        cd = "cd "+ "\""+str(pathlib.Path(self.path2dsmc).absolute())+"\""
        cmd = "@"+self.bash_path + " --login -c " + cd + " " + instr
        f.write_text(cmd)
        return shlex.split(instr)

class SpartaIO:
    ''' Input-Output module for interfacing with SPARTA'''

    SIMID = "simname"
    SIMFILE = "simfile"

    def __init__(self):
        self.writer = Sparta_Writer()
        self.reader = Sparta_Reader()

    def write(self,sim_param,fresh = False):
        ''' Write a simulation (input)'''
        out_params = self.reader.prepare_out(sim_param[SimInput.ID],self.sim_opt)
        return self.writer.write_simulation(sim_param,out_params,fresh)

    def read(self,simid, sim_opt):
        ''' Read a simulation (output)'''
        return self.reader.readsurf(simid, sim_opt)

    def define_simulation(self,sim_opt):
        ''' Define output options applied for all simulations '''
        self.sim_opt = self.reader.set_options(sim_opt)

    def save_simulation(self,name = datetime.today()):
        ''' Save all simulation output files into a folder

        Arguments
        ---------
        name : string - Name of the folder into which to move all the output files 
        '''
        movelist = glob.glob("."+sparta_folder+self.reader.log_name+"_*")
        os.makedirs("."+sparta_folder+name,exist_ok = True)
        for fileid in movelist:
            seg = fileid.split("\\")
            filename = seg[-1]
            seg.insert(-1,name)
            newpath = '\\\\'.join(seg)
            os.rename(fileid,newpath)   

class Sparta_Writer:
    ''' Writes data input for SPARTA exe'''
    base_name = "axi.in"
    example_input = "axi.init"
    prime_input = "axi.init"
    input_path = os.getcwd()+sparta_folder
    token = "§"

    def __init__(self):
        self.clear_folder()
        self.simlist = []

    def write_simulation(self,sim_param,out_params,fresh = False):
        ''' Write simulation inut file

        Arguments
        ---------
        sim_param : dict{SimInput: string} - simulation parameters dictionary
        out_params : dict{SimOutput: string} - simulation output file name dictionary
        --------- Keyword args
        fresh : boolean - Does this require the full simulation file ?
        '''
        self.checkparams(sim_param)
        simid = sim_param.pop(SimInput.ID)
        self.simlist.append(simid)
        simname, siminstr = self.prepare_simfile(simid,fresh)
        with open(self.input_path+simname,"w") as f:
            for k,v in sim_param.items():
                siminstr = siminstr.replace(self.token+k.value,str(v),1)
            for k,v in out_params.items():
                siminstr = siminstr.replace(Sparta_Reader.activation+k.value,"")
                siminstr = siminstr.replace(self.token+k.value,str(v))
            f.write(siminstr)
        sim_param["simname"] = simid
        sim_param["simfile"] = simname
        sim_param.update(out_params)
        return simname

    def checkparams(self,params):
        ''' Check if required input parameters are all defined'''
        missing = []
        for p in SimInput:
            if p not in params:
                missing.append(p)
        if(len(missing) != 0):
            raise Exception("Missing input params, {0} undefined".format(missing))


    def prepare_simfile(self,id,fresh):
        ''' Generate input file according to simulation id'''
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
        ''' Erase all input files from previous simulations'''
        clearlist = glob.glob("."+sparta_folder+self.base_name+"_*")
        clearlist.extend(glob.glob("."+sparta_folder+self.prime_input+"_*"))
        print("Clearing folder for new simulation inputs : ",str(len(clearlist))," fileid(s) removed")
        for clearfile in clearlist:
            os.remove(clearfile)

class Sparta_Reader:
    ''' Reads output files from SPARTA exe'''

    pre = ["v_","c_","f_"]
    log_name = "axi.log"
    surfbase_name = "axi.sout"
    flowbase_name = "axi.fout"
    output_path = os.getcwd()+sparta_folder
    activation = "#§"

    simout_link = {SimOutput.SURF_FILE:surfbase_name,SimOutput.FLOW_FILE:flowbase_name, SimOutput.LOG_FILE:log_name}
    
    def __init__(self):
        self.clear_folder()
        pass
    
    def readsurf(self,simid, sim_opt):
        ''' Read surface data

        Arguments
        ---------
        simid : string - simulation identifier
        sim_opt : dict{SimOptions: string} - simulation output option dictionary

        Returns
        -------
        output dictionary {SimOptions : float[...] array-like}
        '''
        logreader = olog("."+sparta_folder+self.log_name+"_"+simid) # Reads log files only
        outnames = [self.pre[0]+opt.key for opt in sim_opt] # Quantities to be read from log file
        res = logreader.get(*outnames)
        dic = {}
        for i, key in enumerate(sim_opt):
            dic[key] = res[i]
        return dic

    def clear_folder(self):
        ''' Erases all output files from previous simulations '''
        clearlist = glob.glob("."+sparta_folder+self.surfbase_name+"_*")
        clearlist.extend(glob.glob("."+sparta_folder+self.log_name+"_*"))
        clearlist.extend(glob.glob("."+sparta_folder+self.flowbase_name+"_*"))
        print("Clearing folder for new simulation outputs : ",str(len(clearlist))," fileid(s) removed")
        for clearfile in clearlist:
            os.remove(clearfile)

    def set_options(self,opts):
        ''' Processes the options required by the users and returns a full dict of options used for the simulation input files

        Arguments
        ---------
        opts : list[SimOptions] - list of output values that are required to be computed for the user's purpose
        '''
        self.read_opt = opts
        sim_opt = {}
        outdep = {}
        link = {}
        order = 1
        for opt in opts:
            sim_opt[opt] = opt.key
            sim_opt[opt.index] = order
            for dep in opt.dependencies: # Options have dependent enum objects that specify where they will be output into
                if (dep not in outdep) and (dep in self.simout_link): # add a new dependency only if it is implemented and not yet added
                    outdep[dep] = self.simout_link[dep]
                if dep not in self.simout_link: # That means the options must be grouped into the same group-string (for the input file)                    
                    if dep not in link: # Start a new list
                        link[dep] = [opt.key," "]
                    else: # Extend the existing list
                        link[dep].extend([opt.key," "]) 
            order += 1 # The next option is then written in the second place
        sim_opt.update(outdep) # Add the dependencies to the options
        for dep,lopt in link.items():
            sim_opt[dep] = "".join(lopt) # Concatenate the list of string into one single line to be replaced in the input file
        return sim_opt

    def prepare_out(self,simid,sim_opt):
        ''' Prepare the simulation output files'''
        out_params = {}
        out_params.update(sim_opt)
        for opt in sim_opt.keys():
            if opt in self.simout_link:
                out_params[opt] += "_" + simid
        return out_params

    def sort_out(self,prefix = log_name):
        ''' Sort and return the list of output files specified by a prefix name'''
        filelist = glob.glob(prefix+"*")
        self.offset = len(prefix)
        filelist.sort(key=self.find_index)
        return filelist            

    def find_index(self,filename):
        ''' Find the simulation index specifier in a filename string'''
        i = self.offset
        while(not filename[i].isdigit()): i+=1 # Number has more digits
        j = i+1
        while(filename[j].isdigit()): j+=1 # End of the number index
        n = int(filename[i:j]) # From first to last digit
        return n


    def readsurf2(self):
        ''' Unused code'''
        step = [] #this should be the parameter of the object which gives the information of selected timestep
        path = self.output_path
        arraystr = [] #to store information in string
        txt = glob.glob(path)
        #read fileid using dump and separate them by timestep using scatter
        for textfile in txt:
            #find the key name of the fileid such as pressure.surf -> take pressure and put it as a name for scatter arguement
             names = re.split('\\\\',textfile)
             key = re.split('\.',names[-1])
             d = dump(textfile)
             d.scatter(str(key[0]))
        #read the scattered files and prepare for sorting (for .surf files)
        newPath = r'C:\cygwin\home\ASiapan\sparta\bench\pressEng.*'
        txt1 = glob.glob(newPath)
        i = 0
        checkname = 0
        loop = 0
        sortmax = []
        timestep = 100

# Finding max timestep
        for name in txt1:
            st = re.split('\\\\',name)
            fileindex = re.split('\.',st[-1])
            sortmax.append(int(fileindex[-1]))

        sortmax.sort()

#Start collecting data
        while loop == 0:
            for name in txt1:
                st = re.split('\\\\',name)
                fileindex = re.split('\.',st[-1])
        #        Sort fileid according to timestep to collect information in order of timestep in step[]
                if int(fileindex[-1]) == checkname:
                    if int(fileindex[-1]) == max(sortmax):
                        loop =1
                    arrayfl = [] #for float storage
                    arraytmp = []
                    j = 0
                    f = open(name,'r')
                    store = f.read()
                    f.close()
                    arraystr.append(re.split('\n',store))
        #            turn string to float
                    for many in arraystr[i]:
                        tmp = re.split(' ',many)
                        if tmp != ['']:
                            del tmp[2:3]
                            elemno = 0
                            if j >= 9:
                                for elem in tmp:
                                    tmp[elemno] = float(elem)
                                    elemno += 1
                            arraytmp.append(tmp)
                        j += 1
                    del arraytmp[0:9]
        #            sort the surface
                    row = 0
                    lop = 0
                    while lop == 0:
                        if row <= 11:
                            arrayfl.append(arraytmp[row])
                            row += 13
                            arrayfl.append(arraytmp[row])
                            row += 13
                            arrayfl.append(arraytmp[row])
                            row += 12
                            arrayfl.append(arraytmp[row])
                            row -= 37
                            lop = 0
                        else:
                            arrayfl.append(arraytmp[row])
                            row += 13
                            arrayfl.append(arraytmp[row])
                            lop = 1                
                    i += 1
                    checkname += timestep
                    self.step.append(arrayfl)


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