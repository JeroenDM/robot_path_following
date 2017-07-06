
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Joint:
    """ Containts joint information to construct robot """
    def __init__(self, offset, angle, joint_type):
        self.d = offset     #initial value if variable
        self.theta = angle  #initial value if variable
        if joint_type == 'R' or joint_type == 'P':
            self.type = joint_type # 'R' or 'P'
        else:
            raise TypeError("Joint type must be 'R' or 'P'")
        

class Link:
    """ Containts link information to construct robot """
    def __init__(self, length, twist):
        self.a = length
        self.alpha = twist
        
    def __str__(self):
        return "This is a link"

class End_effector:
    """ End effector """
    def __init__(self, position):
        self.p = position # position of end effector in last link frame
        Tee = np.eye(4)
        Tee[0:3, 3] = np.array(position)
        self.T = Tee

class Robot:
    """ Bundle with robot information """
    def __init__(self, links, joints, end_effector):
        self.L = links
        self.J = joints
        self.ee = end_effector
        if len(links) == len(joints):
            self.ndof = len(joints)
        else:
            raise TypeError("Number of links must equal number of joints")
        
    def __str__(self):
        return "This is a robot"
        
    def get_T(self, i, qi):
        """ Get transformation matrix for link
        
        Type must be inferred from imput elements to allow for usage with casadi
        
        Parameters
        ----------
        i : integer
          Index of link for which R is wanted
        q : float or object
          joint value of joint i
          
        returns
        -------
        T : ndarray
          Homogenous transformation matrix of link relative to previous link
        
        """
        from numpy import cos, sin
        
        # link parameters always costant  
        a = self.L[i].a
        alpha = self.L[i].alpha
        
        # joint values, only value of prismatic joint used in this function
        if self.J[i].type == 'R':
          d = self.J[i].d
        else:
          d = qi + self.J[i].d # add initial value
        
        T = np.zeros((4, 4), dtype=type(qi))
        T[-1, -1] = 1
        T[0:3, -1] = np.array([a, -sin(alpha) * d, cos(alpha) * d])
        T[0:3, 0:3] = self.get_R(i, qi)
        
        return T
        
    def get_R(self, i, qi):
        """ Get rotation matrix for a link 
        
        Type must be inferred from imput elements to allow for usage with casadi
        Parameters
        ----------
        i : integer
          Index of link for which R is wanted
        q : float or object
          joint value of joint i
          
        returns
        -------
        R : ndarray
          Rotation matrix of current link relative to previous link
        
        """
        from numpy import cos, sin
        
        # link parameters always costant
        # a = self.L[i].a
        alpha = self.L[i].alpha
        
        # joint parameter, only angle relevant for rotation matrix
        if self.J[i].type == 'R':
          theta = qi + self.J[i].theta # add initial value
        else:
          theta = self.J[i].theta
        
        R = np.zeros((3, 3), dtype=type(theta))
        
        # precompute cos and sin ?
        R[0, ] = np.array([cos(theta), -sin(theta), 0])
        R[1, ] = np.array([sin(theta) * cos(alpha), cos(theta) * cos(alpha), -sin(alpha)])
        R[2, ] = np.array([sin(theta) * sin(alpha), cos(theta) * sin(alpha),  cos(alpha)])
        
        return R
        
    def get_fk_T(self, q):
        """ Get transformation matrix of ee relative to base """
        # figure out a appropriate type for T
        # use type of first joint value
        T = np.eye(4, dtype=type(q[0]))
        for i in range(self.ndof):
            Ti = self.get_T(i, q[i])
            T = T.dot(Ti)
        T = T.dot(self.ee.T)
        return T
        
    def plot2D(self, q):
        """ Draw a statick plot of the robot in 2D """       
        plt.figure()
        plt.axis([-3, 3, -3, 3])
        plt.xlabel('x')
        plt.ylabel('y')
        
        res = self.get_joint_xy(q)
        x = res['x']
        y = res['y']
        plt.plot(x, y, 'ko-')
        
        xee = res['xee']
        yee = res['yee']
        plt.plot(xee, yee, 'ko--')
        
        plt.show()
        
    def get_joint_xy(self, q):
        """ Get the x and y coordinates of joints and end effector """
         # extract x y coordinate from Transform matrices
        x = [0.0]
        y = [0.0]
        T = np.eye(4)
        for i in range(self.ndof):
            Ti = self.get_T(i, q[i])
            T = T.dot(Ti)
            x.append(T[0, 3])
            y.append(T[1, 3])
            
        # add end effector
        xee = x[-1]
        yee = y[-1]
        T = T.dot(self.ee.T)
        xee = np.append(xee, T[0, 3])
        yee = np.append(yee, T[1, 3])
        
        return {'x': x, 'y': y, 'xee': xee, 'yee': yee}
        
        
    def animate(self, joint_values):
        joint_values = np.array(joint_values)
        
        ndof, N = joint_values.shape
        
        if ndof != self.ndof:
            raise ValueError("joint_values must have ndof rows")
        
        l_max = 3.0
        #l_max = sum([ self.L[i].a for i in range(self.ndof) ])
        #l_max += max(self.ee.p)
        
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, autoscale_on=False,
                             xlim=(-l_max, l_max),
                             ylim=(-l_max, l_max))
        ax.grid()
        line, = ax.plot([], [], 'ko-', lw=2)
        line_ee, = ax.plot([], [], 'ko--', lw=2)
        
        def init():
            line.set_data([], [])
            line_ee.set_data([], [])
            return line, line_ee


        def animate(i):            
            res = self.get_joint_xy(joint_values[:, i])            
        
            line.set_data(res['x'], res['y'])
            line_ee.set_data(res['xee'], res['yee'])
            return line, line_ee
        
        ani = animation.FuncAnimation(fig, animate, np.arange(1, N),
                                      interval=50, blit=True, init_func=init)
        
        # ani.save('double_pendulum.mp4', fps=15)
        plt.show()
        
            
            
