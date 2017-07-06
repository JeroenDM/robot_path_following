# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 17:59:09 2017

@author: jeroen
"""

import numpy as np
import casadi as ca
import optistack as op
from Robot import *


#links = []
#links.append(Link(0.0, 0.0)) # base link
#links.append(Link(1.0, 0.0))
#links.append(Link(1.0, 0.0))

#ee = End_effector([1.0, 0.0, 0.0])

#joints = []
#joints.append(Joint(0.0, 0.0, 'R'))
#joints.append(Joint(0.0, 0.0, 'R'))
#joints.append(Joint(0.0, 0.0, 'R'))

#bot = Robot(links, joints, ee)
#print bot
#print bot.get_R(0, 1.5)

#print bot.get_T(0, 1.5)

#print bot.ee.T

#print bot.get_fk_T([0.0, 0.0, 0.0])

#bot.plot2D([1.0, -1.0, 0.2])

#P = 30
#q1 = np.linspace(-1.0, 1.0, P)
#q2 = np.linspace(0.1, -0.1, P)
#q3 = np.linspace(-3.0, 0.0, P)

#jv = np.vstack((q1, q2, q3))

#bot.animate(jv)

# other robot

PI_2 = np.pi / 2

links = []
links.append(Link(0.0, 0.0)) # a0 and alpha0
links.append(Link(1.0, 0.0)) # a1 and alpha1
links.append(Link(0.0, -PI_2)) # a2 amd alpha2
links.append(Link(0.0, PI_2))

ee = End_effector([0.0, 1.0, 0.0])

joints = []
joints.append(Joint(0.0, 0.0, 'R')) # d1 and theta1
joints.append(Joint(0.0, 0.0, 'R'))
joints.append(Joint(1.0, 0.0, 'P')) # d3 and theta3
joints.append(Joint(0.0, 0.0, 'R'))

bot = Robot(links, joints, ee)

bot = Robot(links, joints, ee)
print bot
#print bot.get_R(0, 1.5)

print "Frame link 1 rt 0"
print bot.get_T(0, 0.0)
print "Frame link 2 rt 1"
print bot.get_T(1, 0.0)
print "Frame link 3 rt 2"
print bot.get_T(2, 0.0)

I = np.eye(4)
T1 = bot.get_T(0, 0.0)
T2 = bot.get_T(1, 0.0)
T3 = bot.get_T(2, 0.0)

#print bot.ee.T

print "Frame from ee relative to 0"
print bot.get_fk_T([0.0, 0.0, 0.0, 0.0])

#bot.plot2D([0.0, 0.5, 0.0])

P = 30
q1 = np.linspace(0.0, -0.5, P)
q2 = np.linspace(0.0, 0.0, P)
q3 = np.linspace(-0.5, 0.5, P)
q4 = np.linspace(0.0, -0.6, P)


jv = np.vstack((q1, q2, q3, q4))

bot.animate(jv)

np.array([1.0, 1.0]).norm(
