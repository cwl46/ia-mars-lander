import numpy as np
import matplotlib.pyplot as plt
from math import pow, sqrt

# mass, gravitational constant, initial position and velocity
M = 6.42e23
G = 6.67e-11
r = 3.39e6
x0 = np.array([5e6, 0, 0])
v0_sqr = G * M / np.linalg.norm(x0)
v0_2 = np.array([0, sqrt(v0_sqr), 0])
v0_3 = np.array([0, sqrt(1.7 * v0_sqr), 0])
v0_4 = np.array([0, sqrt(2 * v0_sqr), 0])

# simulation time, timestep and time
t_max = 40000
dt = 1
t_array = np.arange(-t_max, t_max, dt)


def positions(v0):
    # initialise variables for start for Euler Integration
    x_euler = x0
    v_euler = v0

    # initialise empty lists to record trajectories for Euler Integration
    x_list_euler = []
    v_list_euler = []

    for t in t_array:
        # append current state to trajectories
        x_list_euler.append(x_euler)
        v_list_euler.append(v_euler)

        # calculate new position and velocity for Euler Integration
        a_euler = - G * M / pow(np.linalg.norm(x_euler), 3) * x_euler
        x_euler = x_euler + dt * v_euler
        v_euler = v_euler + dt * a_euler

    # initialise states to record trajectories for Verlet Integration
    x_list_verlet = [x_list_euler[0]]
    v_list_verlet = [v_list_euler[0]]
    x_verlet = x_list_euler[1]
    v_verlet = v_list_euler[1]

    for t in t_array[1:]:
        # append current state to trajectories
        x_list_verlet.append(x_verlet)
        v_list_verlet.append(v_verlet)

        # calculate new position and velocity for Verlet Integration
        a_verlet = - G * M / pow(np.linalg.norm(x_verlet), 3) * x_verlet
        x_verlet = 2 * x_verlet - x_list_verlet[-2] + pow(dt, 2) * a_verlet
        v_verlet = (x_verlet - x_list_verlet[-1]) / dt

    # return trajectory lists as arrays
    return np.array(x_list_euler), np.array(x_list_verlet)


# build trajectory arrays
x_array_euler_2, x_array_verlet_2 = positions(v0_2)
x_array_euler_3, x_array_verlet_3 = positions(v0_3)
x_array_euler_4, x_array_verlet_4 = positions(v0_4)


# plot the position-time graph
fig, subplots = plt.subplots(2)
plt.suptitle('Trajectory Plots')
subplots[0].title.set_text('Euler')
subplots[1].title.set_text('Verlet')
subplots[0].plot(x_array_euler_2[:, 0], x_array_euler_2[:, 1],
                 'tab:orange', label='Scenario 2: Circular')
subplots[0].plot(x_array_euler_3[:, 0], x_array_euler_3[:, 1],
                 'tab:green', label='Scenario 3: Elliptical')
subplots[0].plot(x_array_euler_4[:, 0], x_array_euler_4[:, 1],
                 'tab:purple', label='Scenario 4: Hyperbolic')
subplots[1].plot(x_array_verlet_2[:, 0], x_array_verlet_2[:, 1],
                 'tab:orange', label='Scenario 2: Circular')
subplots[1].plot(x_array_verlet_3[:, 0], x_array_verlet_3[:, 1],
                 'tab:green', label='Scenario 3: Elliptical')
subplots[1].plot(x_array_verlet_4[:, 0], x_array_verlet_4[:, 1],
                 'tab:purple', label='Scenario 4: Hyperbolic')

for subplot in subplots:
    subplot.grid()
    subplot.legend()
    subplot.add_patch(plt.Circle((0, 0), r, color='tab:blue'))
plt.show()
