"""
Time optimal path following for a 2 DOF planar robot.
Example from the exercise sessions of Optimization of Mechatronic Systems
"""
import meco_binaries;meco_binaries(cpp_splines='master')
from splines import *
import numpy as np
from casadi import *
import matplotlib.pyplot as plt

""" Robot model definition """
grav = 9.81 # gravitational constant

# robot parameters
m_1  = 1.0  # link mass
m_2  = 1.0
I_1  = 0.5  # link mass moment of inertia
I_2  = 0.5
l_1  = 1.0  # link length
l_2  = 1.0
lc_1 = 0.5  # cog location
lc_2 = 0.5
tau_limit_1 = 5.0  # torque limit
tau_limit_2 = 5.0
q_min       = [-np.pi / 2, -np.pi]
q_max       = [ np.pi / 2,  np.pi]
dq_min      = [-np.pi / 2, -np.pi / 2]
dq_max      = [ np.pi / 2,  np.pi / 2]

# dynamic parameters
mu_11 = m_1
mu_21 = m_1 * lc_1
mu_31 = m_1 * lc_1**2 + I_1

mu_12 = m_2
mu_22 = m_2 * lc_2
mu_32 = m_2 * lc_2**2 + I_2

# robot forwar kinematics
# torque = M * ddq + C * dq + G
def chi(q):
    x = l_1 * cos(q[0]) + l_2 * cos(q[0] + q[1])
    y = l_1 * sin(q[0]) + l_2 * sin(q[0] + q[1])
    return [x, y]
    
# robot dynamics
def M(q):
    return np.array([mu_12 * l_1**2 + mu_22 * 2 * l_1 * cos(q[1]) + mu_31 + mu_32,
                     mu_32 + mu_22 * l_1 * cos(q[1]),
                     mu_32 + mu_22 * l_1 * cos(q[1]),
                     mu_32
                     ]).reshape((2, 2))

def C(q, dq):
    return np.array([-mu_22 * l_1 * sin(q[1]) * dq[1],
                     -mu_22 * l_1 * sin(q[1]) * (dq[0] + dq[1]),
                      mu_22 * l_1 * sin(q[1]) * dq[0],
                      0
                      ]).reshape((2, 2))
        
def G(q):
    return np.array([mu_12 * l_1 * cos(q[0]) + mu_21 * cos(q[0]) + mu_22 * cos(q[0] + q[1]),
                     mu_22 * cos(q[0] + q[1])
                     ]).reshape((2,)) * grav

# path that the robot has to follow, a straight line
# as function of a path variable s [0, 1]
def yp(s):
    start = np.array([1.5, 0]).reshape((2,1))
    end   = np.array([0, 1.5]).reshape((2,1))
    return start * s + end * (1 - s)
    
# define spline basis
degree = 3
knotsint = 10
knots = np.hstack(((degree)*[0],np.linspace(0.,1.,knotsint),(degree)*[1]))
m = BSplineBasis(knots,degree)

""" Optimization variables """
# optimization variables
opti = OptiSpline()

# joint angles
q = [opti.Function(m), opti.Function(m)]
dq = [q[0].derivative(1), q[1].derivative(1)]
ddq = [q[0].derivative(2), q[1].derivative(2)]

T = opti.var() # total motion time

print chi([q[0](0.2), q[1](0.2)])


""" Objective """
obj = T


""" constraints """
con = []

# start and end point in joint space
con.append(q[0](0.0) == 0.0)
con.append(q[0](1.0) == 1.0)
con.append(q[1](0.0) == 0.0)
con.append(q[1](1.0) == 1.5)

con.append(dq[0](0.0) == 0.0)
con.append(dq[0](1.0) == 0.0)
con.append(dq[1](0.0) == 0.0)
con.append(dq[1](1.0) == 0.0)

# joint limits
con.append(q[0] >= q_min[0])
con.append(q[0] <= q_max[0])
con.append(q[1] >= q_min[1])
con.append(q[1] <= q_max[1])

# speed limits
con.append(dq[0] >= dq_min[0] * T)
con.append(dq[0] <= dq_max[0] * T)
con.append(dq[1] >= dq_min[1] * T)
con.append(dq[1] <= dq_max[1] * T)

# torque limits at control points
CP = 22 # numer of control points for constraints
t = np.linspace(0.0, 1.0, CP)
for ti in t:
    q_i = np.array( [q[0](ti), q[1](ti)] )
    dq_i = np.array( [dq[0](ti), dq[1](ti)] )
    ddq_i = np.array( [ddq[0](ti), ddq[1](ti)] )
    tau_i = M(q_i).dot(ddq_i) + C(q_i, dq_i).dot(dq_i) + G(q_i) * np.array([T])**2
    con.append(tau_i[0] <=  tau_limit_1 * T**2)
    con.append(tau_i[0] >= -tau_limit_1 * T**2)
    con.append(tau_i[1] <=  tau_limit_2 * T**2)
    con.append(tau_i[1] >= -tau_limit_2 * T**2)


""" Set initial conditions and solve problem """    
# solve optimisation problem 

sol = opti.solver(obj,con,"ipopt")

q1_init = np.linspace(0.0, 1.0, m.dimension())
q2_init = np.linspace(0.0, 1.5, m.dimension())

sol.value(q[0].coeff(), q1_init)
sol.value(q[1].coeff(), q2_init)

#time = Parameter()
#sol.value(q[0], time)
#sol.value(q[1], time)

sol.solve()

""" Look at solution """
DP = 50 # points to evalutate for plotting
t = np.linspace(0.0, 1.0, DP)
T_ = sol.value(T)
q_ = [sol.value(q[0]), sol.value(q[1])]
q1 = q_[0].list_eval(t)
q2 = q_[1].list_eval(t)
dq1 = q_[0].derivative(1).list_eval(t) / T_
dq2 = q_[1].derivative(1).list_eval(t) / T_
ddq1 = q_[0].derivative(2).list_eval(t) / T_**2
ddq2 = q_[1].derivative(2).list_eval(t) / T_**2


print T_

plt.figure()
plt.title('Joint angles')
plt.plot(t, q1, t, q2)
plt.legend(['q1', 'q2'])
plt.xlabel('Normalized time')
plt.ylabel('Angle [rad]')

plt.figure()
plt.title('Joint speed')
plt.plot(t, dq1, 'r', t, dq2, 'g')
plt.legend(['dq1', 'dq2'])
plt.plot(t, dq_min[0] * np.ones(t.shape), "r--" )
plt.plot(t, dq_max[0] * np.ones(t.shape), "r--" )
plt.plot(t, dq_min[1] * np.ones(t.shape), "g--" )
plt.plot(t, dq_max[1] * np.ones(t.shape), "g--" )
plt.xlabel('Normalized time')
plt.ylabel('Speed [rad/s]')

# calculate torque solution
tau = np.zeros((2, DP))
for i in range(DP):
    q_i   = np.array([q1[i], q2[i]])
    dq_i  = np.array([dq1[i], dq2[i]])
    ddq_i = np.array([ddq1[i], ddq2[i]])
    tau[:, i] = M(q_i).dot(ddq_i) + C(q_i, dq_i).dot(dq_i) + G(q_i)

    
plt.figure()
plt.title('Joint troque')
plt.plot(t, tau[0, :], 'r', t, tau[1, :], 'g')
plt.legend(['tau1', 'tau2'])
plt.plot(t, tau_limit_1 * np.ones(t.shape), "r--" )
plt.plot(t, -tau_limit_1 * np.ones(t.shape), "r--" )
plt.plot(t, tau_limit_2 * np.ones(t.shape), "g--" )
plt.plot(t, -tau_limit_2 * np.ones(t.shape), "g--" )
plt.xlabel('Normalized time')
plt.ylabel('Torque [Nm]')

