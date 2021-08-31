import numpy as np
import matplotlib.pyplot as plt
from math import pow

# mass, gravitational constant, initial position and velocity
M = 6.42e23
G = 6.67e-11
r = 3.39e6
x0 = np.array([3.41e6, 0, 0])
v0 = np.array([0, 0, 0])

# simulation time, timestep and time
t_max = 100
dt = 0.1
t_array = np.arange(0, t_max, dt)

# initialise variables for start
x_euler = x0
v_euler = v0

# initialise empty lists to record trajectories for Euler Integration
x_list_euler = []
v_list_euler = []

for t in t_array:
    # append current state to trajectories
    x_list_euler.append(x_euler)
    v_list_euler.append(v_euler)

    # calculate new position and velocity - Euler Integration
    a_euler = - G * M / pow(np.linalg.norm(x_euler), 3) * x_euler
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
    a_verlet = - G * M / pow(np.linalg.norm(x_verlet), 3) * x_verlet
    x_verlet = 2 * x_verlet - x_list_verlet[-2] + pow(dt, 2) * a_verlet
    v_verlet = (x_verlet - x_list_verlet[-1]) / dt


# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2)
x_array_euler = np.array(
    [(np.linalg.norm(x_euler) - r) for x_euler in x_list_euler])
v_array_euler = np.array(
    [np.linalg.norm(v_euler) for v_euler in v_list_euler])
x_array_verlet = np.array(
    [(np.linalg.norm(x_verlet) - r) for x_verlet in x_list_verlet])
v_array_verlet = np.array(
    [np.linalg.norm(v_verlet) for v_verlet in v_list_verlet])

# plot the position-time graph
fig, subplots = plt.subplots(2, 2)
plt.suptitle('Scenario 1: Vertical descent')
subplots[0][0].plot(t_array, x_array_euler)
subplots[0][0].title.set_text('Euler: Height (m)')
subplots[0][1].plot(t_array, v_array_euler, 'tab:orange')
subplots[0][1].title.set_text('Euler: Velocity (m/s)')
subplots[1][0].plot(t_array, x_array_verlet, 'g')
subplots[1][0].title.set_text('Verlet: Height (m)')
subplots[1][1].plot(t_array, v_array_verlet, 'tab:purple')
subplots[1][1].title.set_text('Verlet: Velocity (m/s)')
for subplotrow in subplots:
    for subplot in subplotrow:
        subplot.grid()
plt.show()
