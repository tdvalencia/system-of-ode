###########################################
# Solves ODE equations using methods from
# Differential Equations course (CALC 4)
###########################################

import math
import numpy as np

def dydt(t, y):
    return 1.1 - y + y**3

def euler_method(y0, t0, tn, h):
    t_values = [t0]
    y_values = [y0]
    while t_values[-1] < tn:
        y_next = y_values[-1] + h * dydt(t_values[-1], y_values[-1])
        y_values.append(y_next)
        t_values.append(t_values[-1] + h)
    return t_values, y_values

def improved_euler(y0, t0, tn, h):
    t_values = [t0]
    y_values = [y0]
    while t_values[-1] < tn:
        y_next = y_values[-1] + (h/2) * ( dydt(t_values[-1], y_values[-1])
                                            + dydt(t_values[-1] + h, y_values[-1]
                                            + h*dydt(t_values[-1], y_values[-1])) )
        y_values.append(y_next)
        t_values.append(t_values[-1] + h)
    return t_values, y_values

def runge_kutta_4th(y0, t0, tn, h):
    t_values = [t0]
    y_values = [y0]
    while t_values[-1] < tn:
        k1 = h * dydt(t_values[-1], y_values[-1])
        k2 = h * dydt(t_values[-1] + (h/2), y_values[-1] + (k1/2))
        k3 = h * dydt(t_values[-1] + (h/2), y_values[-1] + (k2/2))
        k4 = h * dydt(t_values[-1] + h, y_values[-1] + k3)

        y_next = y_values[-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

        y_values.append(y_next)
        t_values.append(t_values[-1] + h)
    return t_values, y_values

def improved_euler_tolerance(y0, x0, c, tol, M):
    z = y0
    flag = False
    for m in range(M):
        N = 2**m

        ## Improved Euler method ##
        h = (c-x0)/N; x = x0; y = y0
        for i in range(N):
            F = dydt(x,y)
            G = dydt(x+h,y+h*F)

            x = x+h
            y = y+h*(F+G)/2
        ###########################
        print(f"h={h:.8f}, y={y:.8f}, error={y-z}")
        if abs(y-z) < tol:
            flag = True
            break
        z = y
    if flag:
        print(f"phi(c) ~= {y} with tolerance {tol}")
    else:
        print(f"phi(c) ~= {y} but may not be in tolerance {tol}")

if __name__ == '__main__':
    # Initial conditions
    y0 = 0
    t0 = 0
    tn = 0.9

    # Applying Euler's or Runge-Kutta method
    # h = 0.02
    # t_values, y_values = runge_kutta_4th(y0, t0, tn, h)

    # Output the results
    # for i in range(len(t_values)):
    #     print(f"t = {t_values[i]:.3f}, y = {y_values[i]:.7f}")

    # Applying improved Euler's method with tolerance
    tol = 0.007
    M = 100
    improved_euler_tolerance(y0, t0, tn, tol, M)

