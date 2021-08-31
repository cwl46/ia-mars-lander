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
x_verlet = x0
x_verlet_prev = x0
v_verlet = v0

# initialise empty lists to record trajectories
x_list_euler = []
v_list_euler = []
x_list_verlet = []
v_list_verlet = []

for t in t_array:
    # append current state to trajectories
    x_list_euler.append(x_euler)
    v_list_euler.append(v_euler)
    x_list_verlet.append(x_verlet)
    v_list_verlet.append(v_verlet)

    # calculate new position and velocity
    a_euler = - G * M / pow(np.linalg.norm(x_euler), 3) * x_euler
    x_euler = x_euler + dt * v_euler
    v_euler = v_euler + dt * a_euler

    a_verlet = - G * M / pow(np.linalg.norm(x_verlet), 3) * x_verlet
    x_verlet = 2 * x_verlet - x_verlet_prev + pow(dt, 2) * a_verlet
    x_verlet_prev = x_list_verlet[-1]
    v_verlet = (x_verlet - x_list_verlet[-1]) / dt

# convert trajectory lists into arrays
h_array_euler = np.linalg.norm(np.array(x_list_euler), axis=1) - r
v_array_euler = np.linalg.norm(np.array(v_list_euler), axis=1)
h_array_verlet = np.linalg.norm(np.array(x_list_verlet), axis=1) - r
v_array_verlet = np.linalg.norm(np.array(v_list_verlet), axis=1)

# plot the position-time graph
fig, subplots = plt.subplots(2, 2)
plt.suptitle('Scenario 1: Vertical descent')
subplots[0, 0].plot(t_array, h_array_euler, 'tab:blue')
subplots[0, 0].title.set_text('Euler: Height (m)')
subplots[0, 1].plot(t_array, v_array_euler, 'tab:orange')
subplots[0, 1].title.set_text('Euler: Velocity (m/s)')
subplots[1, 0].plot(t_array, h_array_verlet, 'tab:green')
subplots[1, 0].title.set_text('Verlet: Height (m)')
subplots[1, 1].plot(t_array, v_array_verlet, 'tab:purple')
subplots[1, 1].title.set_text('Verlet: Velocity (m/s)')
for subplotrow in subplots:
    for subplot in subplotrow:
        subplot.grid()
plt.show()
