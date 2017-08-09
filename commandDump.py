# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 12:15:49 2017

@author: ninja_000
"""
import os
import re
from dump import dump
import glob
import matplotlib.pyplot as mpl
import numpy as np

class Sparta_Reader:
    
    def __init__(self):
        
    
    def readsurf(self):
        step = [] #this should be the parameter of the object which gives the information of selected timestep
        path = r'D:\Program\cygwin\home\ninja_000\sparta\bench\*.surf'
        arraystr = [] #to store information in string
        txt = glob.glob(path)
        #read file using dump and separate them by timestep using scatter
        for textfile in txt:
            #find the key name of the file such as pressure.surf -> take pressure and put it as a name for scatter arguement
             names = re.split('\\\\',textfile)
             key = re.split('\.',names[-1])
             d = dump(textfile)
             d.scatter(str(key[0]))
        #read the scattered files and prepare for sorting (for .surf files)
        newPath = r'D:\Program\cygwin\home\ninja_000\sparta\tools\pizza\pressEng.*'
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
        #        Sort file according to timestep to collect information in order of timestep in step[]
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
#
#arr = np.array(self.step)/2.9218564839510575
#mpl.plot(arr[10,:,0])
#
#
#pathgrid = r'D:\Program\cygwin\home\ninja_000\sparta\bench\*.grid'
#txtg = glob.glob(pathgrid)
#for textfile in txtg:
#  names = re.split('\\\\',textfile)
#  key = re.split('\.',names[-1])
#  dum = dump(textfile)
#  dum.scatter(str(key[0]))