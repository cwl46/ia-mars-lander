import numpy as np
import matplotlib.pyplot as plt
euler_results = np.loadtxt('trajectories_euler.txt')
verlet_results = np.loadtxt('trajectories_verlet.txt')

fig, subplots = plt.subplots(2, sharex=True)
subplots[0].plot(euler_results[:, 0], euler_results[:, 1], label='x (m)')
subplots[0].plot(euler_results[:, 0], euler_results[:, 2], label='v (m/s)')
subplots[0].title.set_text('Euler Integration')
subplots[1].plot(verlet_results[:, 0], verlet_results[:, 1], label='x (m)')
subplots[1].plot(verlet_results[:, 0], verlet_results[:, 2], label='v (m/s)')
subplots[1].title.set_text('Verlet Integration')
for subplot in subplots:
    subplot.grid()
    subplot.legend()
plt.show()
