import numpy as np
from matplotlib import pyplot as plt

# Function Definitions
def dxdt(x, y):
    return x + y

def dydt(x, y):
    return x - y

# Euler's method: deltay ~= (dy/dx)deltax
def eulerdt(x0, y0, t0, tf, deltat):
    Lx = []; Ly = []
    x = x0; y = y0
    for t in np.arange(t0, tf, deltat):
        Lx.append(x); Ly.append(y)

        # Recalculate x and y
        x = x + dxdt(x,y) * deltat
        y = y + dydt(x,y) * deltat

    return Lx, Ly


# Run and Plot
if __name__ == '__main__':

    # Initial Values
    deltat = 0.0001 # step size
    x0 = -1; y0 = 1.5
    t0 = 0; tf = 1.5

    # Euler Method
    Lx, Ly = eulerdt(x0, y0, t0, tf, deltat)
    # print(f'Lx: {Lx}'); print(f'Ly: {Ly}')
    plt.plot(Lx, Ly)

    # Slope Field
    nx, ny = 0.2, 0.2
    x = np.arange(-4, 4, nx)
    y = np.arange(-3, 3, ny)

    X, Y = np.meshgrid(x, y)

    dy = dydt(X,Y)
    dx = dxdt(X,Y)

    plt.quiver(X, Y, dx, dy, color='purple')

    # Plot Options
    plt.grid()
    plt.xlim([-3,3])
    plt.ylim([-2,2])
    plt.title('y(x)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
