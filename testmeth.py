# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 10:38:59 2017

@author: ASiapan

Tester class for various methods of the code
"""

import unittest
import numpy as np
import methutil as methu
from mathutils import Quaternion

class TestMethu(unittest.TestCase):

    def setUp(self):
        pass

    def test_anglevec(self):
        vecref = np.array([0,0,1])
        quat = Quaternion(vecref,0.1)
        vec1 = np.array([1,0,0])
        vec2 = np.array([0,1,0])
        ang = methu.angle_vec(vec1,vec2,vecref)[0]
        revang = methu.angle_vec(vec2,vec1,vecref)[0]
        print('{0} x {1} -> {2}'.format(vec1,vec2,methu.angle_vec(vec1,vec2,vecref)[0]))
        self.assertEqual(ang,np.pi/2,msg=ang)
        self.assertEqual(revang,3*np.pi/2,msg=revang)

        quat.angle = np.pi/2
        vec2 = np.array(quat.to_matrix()).dot(vec2)
        print('{0} x {1} -> {2}'.format(vec1,vec2,methu.angle_vec(vec1,vec2,vecref)[0]))
#        self.assertEqual(methu.angle_vec(vec1,vec2,vecref)[0],np.pi,msg=methu.angle_vec(vec1,vec2,vecref)[0])

        quat.angle = np.pi/2
        vec2 = np.array(quat.to_matrix()).dot(vec2)
        print('{0} x {1} -> {2}'.format(vec1,vec2,methu.angle_vec(vec1,vec2,vecref)[0]))
#        self.assertEqual(methu.angle_vec(vec1,vec2,vecref)[0],np.pi*3/2,msg=methu.angle_vec(vec1,vec2,vecref)[0])

unittest.main()