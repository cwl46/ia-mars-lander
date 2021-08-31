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
t_max = 1500
dt = 1.9999
t_array = np.arange(0, t_max, dt)

# initialise variables for start
x_euler = x0
v_euler = v0

# initialise empty lists to record trajectories
x_list_analytical = []
v_list_analytical = []
x_list_euler = []
v_list_euler = []

for t in t_array:
    # append current state to trajectories
    x_list_analytical.append(x_analytical(omega, phi, A, t))
    v_list_analytical.append(v_analytical(omega, phi, A, t))
    x_list_euler.append(x_euler)
    v_list_euler.append(v_euler)

    # calculate new position and velocity - Euler Integration
    a_euler = - k * x_euler / m
    x_euler = x_euler + dt * v_euler
    v_euler = v_euler + dt * a_euler

# initialise initial states to record trajectories for Verlet Integration
x_list_verlet = [x_list_euler[0]]
v_list_verlet = [v_list_euler[0]]
x_verlet = x_list_euler[1]
v_verlet = v_list_euler[1]

for t in t_array[1:]:
    # append current state to trajectories
    x_list_verlet.append(x_verlet)
    v_list_verlet.append(v_verlet)

    # calculate new position and velocity - Verlet Integration
    a_verlet = - k * x_verlet / m
    x_verlet = 2 * x_verlet - x_list_verlet[-2] + math.pow(dt, 2) * a_verlet
    v_verlet = (x_verlet - x_list_verlet[-1]) / dt


# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2)
x_array_analytical = np.array(x_list_analytical)
v_array_analytical = np.array(v_list_analytical)
x_array_euler = np.array(x_list_euler)
v_array_euler = np.array(v_list_euler)
x_array_verlet = np.array(x_list_verlet)
v_array_verlet = np.array(v_list_verlet)

# plot the position-time graph
fig, subplots = plt.subplots(3, sharex=True)
subplots[0].plot(t_array, x_array_analytical, label='x (m)')
subplots[0].plot(t_array, v_array_analytical, label='v (m/s)')
subplots[0].title.set_text('Analytical Solution dt={}'.format(dt))
subplots[1].plot(t_array, x_array_euler, label='x (m)')
subplots[1].plot(t_array, v_array_euler, label='v (m/s)')
subplots[1].title.set_text('Euler Integration dt={}'.format(dt))
subplots[2].plot(t_array, x_array_verlet, label='x (m)')
subplots[2].plot(t_array, v_array_verlet, label='v (m/s)')
subplots[2].title.set_text('Verlet Integration dt={}'.format(dt))
for subplot in subplots:
    subplot.grid()
    subplot.legend()
plt.show()
