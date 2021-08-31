#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main()
{

  // declare variables
  double m, k, x0, v0, x_euler, v_euler, a_euler, x_verlet, v_verlet, a_verlet, t_max, dt, t;
  vector<double> t_list, x_euler_list, v_euler_list, x_verlet_list, v_verlet_list;

  // mass, spring constant, initial position and velocity
  m = 1;
  k = 1;
  x0 = 0;
  v0 = 1;

  // simulation time and timestep
  t_max = 100;
  dt = 0.1;

  // initialise variables for start for Euler Integration
  x_euler = x0;
  v_euler = v0;

  // Euler integration
  for (t = 0; t <= t_max; t = t + dt)
  {
    // append current state to trajectories
    t_list.push_back(t);
    x_euler_list.push_back(x_euler);
    v_euler_list.push_back(v_euler);

    // calculate new position and velocity
    a_euler = -k * x_euler / m;
    x_euler = x_euler + dt * v_euler;
    v_euler = v_euler + dt * a_euler;
  }

  // initialise states for start for Verlet Integration
  x_verlet_list.push_back(x_euler_list[0]);
  v_verlet_list.push_back(v_euler_list[0]);
  x_verlet = x_euler_list[1];
  v_verlet = v_euler_list[1];

  for (t = dt; t <= t_max; t = t + dt)
  {
    // append current state to trajectories
    x_verlet_list.push_back(x_verlet);
    v_verlet_list.push_back(v_verlet);

    // calculate new position and velocity
    a_verlet = -k * x_verlet / m;
    x_verlet = 2 * x_verlet - x_verlet_list[x_verlet_list.size() - 2] + pow(dt, 2) * a_verlet;
    v_verlet = (x_verlet - x_verlet_list[x_verlet_list.size() - 1]) / dt;
  }

  // Write the trajectories to file
  ofstream fout;
  fout.open("trajectories_euler.txt");
  if (fout)
  { // file opened successfully
    for (int i = 0; i < t_list.size(); i = i + 1)
    {
      fout << t_list[i] << ' ' << x_euler_list[i] << ' ' << v_euler_list[i] << endl;
    }
  }
  else
  { // file did not open successfully
    cout << "Could not open trajectory file for writing" << endl;
  }
  fout.close();
  fout.open("trajectories_verlet.txt");
  if (fout)
  { // file opened successfully
    for (int i = 0; i < t_list.size(); i = i + 1)
    {
      fout << t_list[i] << ' ' << x_verlet_list[i] << ' ' << v_verlet_list[i] << endl;
    }
  }
  else
  { // file did not open successfully
    cout << "Could not open trajectory file for writing" << endl;
  }

  /* The file can be loaded and visualised in Python as follows:

  import numpy as np
  import matplotlib.pyplot as plt
  results = np.loadtxt('trajectories.txt')
  plt.figure(1)
  plt.clf()
  plt.xlabel('time (s)')
  plt.grid()
  plt.plot(results[:, 0], results[:, 1], label='x (m)')
  plt.plot(results[:, 0], results[:, 2], label='v (m/s)')
  plt.legend()
  plt.show()

  */
}
