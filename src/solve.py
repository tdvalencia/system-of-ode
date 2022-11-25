import numpy as np
from matplotlib import pyplot as plt

STEP_SIZE = deltat = 0.5
Lx = []; Ly = []
x0 = 0; y0 = 3

def dxdt(x, y):
    return x + y

def dydt(x, y):
    return x - y

def eulerdt(t0, tf):
    x = x0; y = y0
    for t in np.arange(t0, tf, deltat):
        Lx.append(x); Ly.append(y)
        x = x + dxdt(x, y) * deltat
        y = y + dydt(x, y) * deltat

if __name__ == '__main__':
    eulerdt(-1, 1)
    plt.plot(Lx, Ly)
    plt.show()