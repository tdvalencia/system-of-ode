###########################
# Runge-Kutta 4th order method given
# 2nd order derivative
###########################

import numpy as np
import math

# g(t,y,Dy) is equivalent to y'' where:
#   - t is independent variable
#   - y is dependent
#   - Dy is derivative of y (y')
def g(t,y,Dy):
    return -2*np.sin(y)

# f(t, y, Dy) is equivalent to y' where:
#   - t is independent var
#   - y is dependent
#   - Dy is derivative of y (y')
def f(t,y,Dy):
    return Dy

def vectorize_runge_kutta(t0,tn,y0,Dy0,h):
    t_values = [t0]
    y_values = [y0]
    Dy_values = [Dy0]

    while t_values[-1] < tn:

        # Calculate values for k1;1 and k2;1
        k11 = h * f(t_values[-1], y_values[-1], Dy_values[-1])
        k21 = h * g(t_values[-1], y_values[-1], Dy_values[-1])

        # k1;2 and k2;2
        k12 = h * f(t_values[-1] + (h/2), y_values[-1] + (k11/2),
                      Dy_values[-1] + (k21/2))
        k22 = h * g(t_values[-1] + (h/2), y_values[-1] + (k11/2),
                      Dy_values[-1] + (k21/2))

        # k1;3 and k2;3
        k13 = h * f(t_values[-1] + (h/2), y_values[-1] + (k12/2),
                      Dy_values[-1] + (k22/2))
        k23 = h * g(t_values[-1] + (h/2), y_values[-1] + (k12/2),
                      Dy_values[-1] + (k22/2))

        # k1;4 and k2;4
        k14 = h * f(t_values[-1] + h, y_values[-1] + k13,
                    Dy_values[-1] + k23)
        k24 = h * g(t_values[-1] + h, y_values[-1] + k13,
                    Dy_values[-1] + k23)

        # Combining k's for each Runge-Kutta Calculation
        y_next = y_values[-1] + (1/6)*(k11 + 2*k12 + 2*k13 + k14)
        Dy_next = Dy_values[-1] + (1/6)*(k21 + 2*k22 + 2*k23 + k24)

        y_values.append(y_next)
        Dy_values.append(Dy_next)
        t_values.append(t_values[-1] + h)
    return t_values, y_values, Dy_values

if __name__ == '__main__':
    # Initial conditions
    t0 = 0
    tn = 3
    y0 = 1
    Dy0 = 0

    # Applying Euler's or Runge-Kutta method
    h = 0.02
    t_values, y_values, Dy_values = vectorize_runge_kutta(t0,tn,y0,Dy0,h)

    # Output the results
    for i in range(len(t_values)):
        print(f"t = {t_values[i]:.3f}, y = {y_values[i]:.7f}, Dy = {Dy_values[i]:.7f}")
