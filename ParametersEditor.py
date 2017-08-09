# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 13:54:15 2017

@author: ninja_000
"""
def ParamInput(nv,T,Tvib,v,path):
#Locating Data File
#path = r'C:\Users\ninja_000\Desktop\DS2VD.dat'

#Retrieve old data
    file = open(path,'r')

    oldParam = file.read()

    file.close

#editing data in old data
    arrParam = oldParam.splitlines(500)
#    print(arrParam[399:403])
    densityNew = nv #input('Density (in format of x.xxxxxxe+xx) = ')
    arrParam[399] = str(densityNew) + '\n'
    tempNew = T #input('Temperature (K) = ')
    arrParam[400] = str(tempNew) + '\n'
    vibtempNew = Tvib #input('Temp (vib)(K) = ')
    arrParam[401] = str(vibtempNew) + '\n'
    xveloNew = v #input('Velocity Component in x-direction(m/s) = ')
    arrParam[402] = str(xveloNew) + '\n'

#recompile the string
    newParam = ''
    for index in range(len(arrParam)):
        newParam = newParam + arrParam[index]

#    print(newParam)

#Rewrite the file with new input
    file = open(path,'w')
    file.write(newParam)
    file.close()
