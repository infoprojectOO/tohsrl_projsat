# -*- coding: utf-8 -*-
"""
Created on Thu May 11 13:31:42 2017

@author: ASiapan
"""
import math as m

global meanatm 
meanatm = { 
100 : (184.0160,	5.08E-07),
120	: (374.9715,	1.80E-08),
140	: (635.5703,	3.26E-09),
160	: (787.5532,    1.18E-09),
180	: (877.6729,	5.51E-10),
200	: (931.2806,	2.91E-10),
220	: (963.2701,	1.66E-10),
240	: (982.4191,	9.91E-11),
260	: (993.9173,	6.16E-11),
280	: (1000.8427,	3.94E-11),
300	: (1005.0267,	2.58E-11),
320	: (1007.5620,	1.72E-11),
340	: (1009.1030,	1.16E-11),
360	: (1010.0423,	7.99E-12),
380	: (1010.6166,	5.55E-12),
400	: (1010.9688,	3.89E-12),
420 : (1011.1853, 	2.75E-12)
}

def densityold(alt):
    alt_ref = 71.0 #km
    dst_ref = 6.4*10**(-5) #kg/m^3
    T_ref = 214.65 #K
    M = 0.0289644 #kg/mol
    R = 8.3144 #J/K
    g0 = 9.81 #m/s^2
    
    return dst_ref*math.e**(-(alt-alt_ref)*1000.0*g0*M/(T_ref*R))

def density(alt):
	if(alt<=0):
		#return 10**6
		raise ValueError('Negative Altitude')
	i_measure = 1
	row = int(alt/20) - 5

	h1 = (row+5)*20
	if h1<100:
		h1=100
	elif h1>400:
		h1 = 400
	h2 = h1+20
	if(h2>420):
		h2=420
	d1 = meanatm[h1][i_measure]
	d2 = meanatm[h2][i_measure]

	B = (m.log(d2)-m.log(d1))*0.05
	A = d1/m.e**(B*h1)

	return A*m.e**(B*(alt))
