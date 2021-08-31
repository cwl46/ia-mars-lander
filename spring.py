import numpy as np
import matplotlib.pyplot as plt
import math


def analytical_params(m, k, x0, v0):
    omega = math.sqrt(k / m)
    try:
        phi = - math.atan(v0 / x0 / omega)
        A = x0 / math.cos(phi)
    except ZeroDivisionError:
        phi = math.copysign(math.pi/2, - v0 / omega)
        A = - v0 / omega / math.sin(phi)
    return omega, phi, A


def x_analytical(omega, phi, A, t):
    return A * math.cos(omega * t + phi)


def v_analytical(omega, phi, A, t):
    return - A * omega * math.sin(omega * t + phi)


# mass, spring constant, initial position and velocity
m = 1
k = 1
x0 = 0
v0 = 1
omega, phi, A = analytical_params(m, k, x0, v0)

# simulation time, timestep and time
t_max = 100
dt = 0.1
t_array = np.arange(0, t_max, dt)

# initialise variables for start
x_euler = x0
v_euler = v0
x_verlot = x0
v_verlot = v0
x_verlot_prev = x0
v_verlot_prev = v0

# initialise empty lists to record trajectories
x_list_analytical = []
v_list_analytical = []
x_list_euler = []
v_list_euler = []
x_list_verlot = []
v_list_verlot = []

for t in t_array:
    # calculate new position and velocity - Analytical Soluiton

    # append current state to trajectories
    x_list_analytical.append(x_analytical(omega, phi, A, t))
    v_list_analytical.append(v_analytical(omega, phi, A, t))
    x_list_euler.append(x_euler)
    v_list_euler.append(v_euler)
    x_list_verlot.append(x_verlot)
    v_list_verlot.append(v_verlot)

    # calculate new position and velocity - Euler Integration
    a = - k * x_euler / m
    x_euler = x_euler + dt * v_euler
    v_euler = v_euler + dt * a

# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2)
x_array_analytical = np.array(x_list_analytical)
v_array_analytical = np.array(v_list_analytical)
x_array_euler = np.array(x_list_euler)
v_array_euler = np.array(v_list_euler)

# plot the position-time graph
plt.figure(1)
plt.clf()
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, x_array_euler, label='x (m)')
plt.plot(t_array, v_array_euler, label='v (m/s)')
plt.legend()
plt.show()
