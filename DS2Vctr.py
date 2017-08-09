# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:00:31 2017

@author: ninja_000
"""

import subprocess

#path = 'D:\Program\WinPython-64bit-3.5.3.1Qt5\WinPython Command Prompt'
#subprocess.check_output([path],shell=True)
proc = subprocess.Popen(['D:\Program\WinPython-64bit-3.5.3.1Qt5\WinPython Command Prompt'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = proc.communicate(timeout=1)
#i = 0
#for x in i:
#    
#stdout
#proc.wait()
#print( proc.returncode)
#print('done')