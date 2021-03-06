# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from util import *



q0 = np.array([0.1, -0.3]).reshape(2, 1)
x0 = fk(q0)
J0 = J(q0)

Np = 10
xp = np.vstack((np.linspace(-1.5, 1.5, Np),
                np.ones(Np)
              ))

plt.figure()
plt.plot(xp[0, :], xp[1, :], 'go')
plt.show()
