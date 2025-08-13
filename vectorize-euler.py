######################
# Euler method given
# 2nd order derivative
######################

import numpy as np
import math

# Initial conditions
t0 = 0
tn = 2
y0 = 2
Dy0 = 0
h = 0.5

# g(t,y,Dy) is equivalent to y'' where:
#   - t is independent variable
#   - y is dependent
#   - Dy is derivative of y (y')
def g(t,y,Dy):
    return -2*t*Dy - 4*y

def vectorize_euler(t0,tn,y0,Dy0,h):
    t_values = [t0]
    y_values = [y0]
    Dy_values = [Dy0]

    while t_values[-1] < tn:
        y_next = y_values[-1] + h*Dy_values[-1]
        Dy_next = Dy_values[-1] + h*g(t_values[-1], y_values[-1],
                                      Dy_values[-1])

        y_values.append(y_next)
        Dy_values.append(Dy_next)
        t_values.append(t_values[-1] + h)
    return t_values, y_values, Dy_values

if __name__ == '__main__':
    # Applying Euler's method
    t_values, y_values, Dy_values = vectorize_euler(t0,tn,y0,Dy0,h)

    # Output the results
    for i in range(len(t_values)):
        print(f"t = {t_values[i]:.3f}, y = {y_values[i]:.7f}, Dy = {Dy_values[i]:.7f}")
