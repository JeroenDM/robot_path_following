from numpy.linalg import norm, pinv
from numpy import cos, sin, array, pi, dot, cross, arctan2
from scipy.optimize import root
import matplotlib.pyplot as plt

l1 = 1.0
l2 = 1.0

def fk(q, all_links = False):
  
    x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1])
    y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1])
    if all_links:
        x1 = l1 * cos(q[0])
        y1 = l1 * sin(q[0])
        return array([x, y]), array([x1, y1])
    else:
        return array([x, y])
  
def J(q):
  from numpy import cos, sin, array
  j11 = -l1 * sin(q[0]) - l2 * sin(q[0] + q[1])
  j12 = -l2 * sin(q[0] + q[1])
  j21 =  l2 * cos(q[0] + q[1])
  j22 =  l1 * cos(q[0]) + l2 * cos(q[0] + q[1])
  return array([j11, j12, j21, j22]).reshape(2, 2)
  
def rrplot(q):
  from numpy import cos, sin
  x1 = l1 * cos(q[0])
  y1 = l1 * sin(q[0])
  x2 = l1 * cos(q[0]) + l2 * cos(q[0] + q[1])
  y2 = l1 * sin(q[0]) + l2 * sin(q[0] + q[1])
  x = [ 0.0, x1[0], x2[0] ]
  y = [ 0.0, y1[0], y2[0] ]
  
  plt.figure()
  plt.axis([-3, 3, -3, 3])
  plt.plot(x, y, 'ko-')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.show()
    
    
def dist(p1, p2, a):
    d = norm(p1 - p2)
    if d <= 1e-12:
        raise ValueError("Line must be given by two different points, d= " + str(d))
    A = abs( (p2[1] - p1[1]) * a[0] - (p2[0] - p1[0]) * a[1] + p2[0] * p1[1] - p2[1] * p1[0] )
    return A / d

def mid_point(qv, xv, i):
    q_mid = (qv[i] + qv[i+1]) / 2
    x_mid = fk(q_mid)
    print i
    print xv
    e = dist(xv[i], xv[i+1], x_mid)
    return q_mid, x_mid, e

# add point in joint and path points array
# WARNING: xsol and qsol are modified by the function
def add_point(xsol, qsol, xin, qin, ind):
    xsol.insert(ind, xin)
    qsol.insert(ind, qin)
    
def refine_grid(xsol, qsol, i):
    x_new = (xsol[i] + xsol[i+1]) / 2
    sol = root(lambda q : x_new - fk(q), qsol[i])
    q_new = sol['x']
    
    return x_new, q_new

def angle(a, b, na, nb):
    """ At this point it is checked that the norm is not close to zeros"""
    cos_angle = dot(a, b) / (na * nb)
    sin_angle = cross(a, b) / (na * nb)
    
    return arctan2(sin_angle, cos_angle)

def mass_center(q1, q2, q3, w):
    return (q1 + w * q2 + q3) / (2.0 + w)
  

def newton(x, q0, alpha = 0.1, max_it = 1000, tol = 1e-6):
  it_counter = max_it
  
  qk = q0
  xk = fk(q0)
  Jk = J(q0)
  
  while(it_counter > 0):
    it_counter -= 1
    
    Jk_inv = pinv(Jk)
    q_new = qk + alpha * Jk_inv.dot(x - xk)
    
    xk = fk(q_new)
    if norm(xk - x) <= tol:
      return {'conv': True, 'q': q_new}
    
    qk = q_new
    Jk = J(qk)
    
  return {'conv': False}

def jac_transpose(x, q0, alpha = 0.1, max_it = 1000, tol = 1e-6):
  it_counter = max_it
  
  qk = q0
  xk = fk(q0)
  Jk = J(q0)
  
  while(it_counter > 0):
    it_counter -= 1
    
    Jk_t = Jk.T
    q_new = qk + alpha * Jk_t.dot(x - xk)
    
    xk = fk(q_new)
    if norm(xk - x) <= tol:
      return {'conv': True, 'q': q_new}
    
    qk = q_new
    Jk = J(qk)
    
  return {'conv': False}
