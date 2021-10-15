// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"
#include "math.h"

void autopilot(void)
// Autopilot to adjust the engine throttle, parachute and attitude control
{
  double K_h = 0.018;
  double K_p = 0.9;
  double delta = 0.6;

  double descent_rate = velocity * position.norm();
  vector3d radial_velocity = descent_rate * position.norm();
  vector3d tangential_velocity = velocity - radial_velocity;
  double altitude = position.abs() - MARS_RADIUS;
  double ground_speed = sqrt(velocity.abs2() - pow(descent_rate, 2));

  if (ground_speed > 1.0)
  {
    attitude_stabilization(-tangential_velocity.norm());
    throttle = K_p * (ground_speed - 1.0);
  }
  else
  {
    stabilized_attitude = true;
    attitude_stabilization(position.norm());

    double descent_rate = velocity * position.norm();
    double altitude = position.abs() - MARS_RADIUS;
    double error_term = -(0.5 + K_h * altitude + descent_rate);
    throttle = delta + K_p * error_term;
  }
}

vector3d drag_without_parachute(vector3d position)
{
  double lander_area = M_PI * LANDER_SIZE * LANDER_SIZE;
  return -0.5 * atmospheric_density(position) * DRAG_COEF_LANDER * lander_area * velocity.abs2() * velocity.norm();
}

vector3d drag_with_parachute(vector3d position)
{
  double parachute_area = 5.0 * 2.0 * LANDER_SIZE * 2.0 * LANDER_SIZE;
  return -0.5 * atmospheric_density(position) * DRAG_COEF_CHUTE * parachute_area * velocity.abs2() * velocity.norm();
}

vector3d euler_update(vector3d net_acceleration)
{
  vector3d new_position = position + delta_t * velocity;
  velocity = velocity + delta_t * net_acceleration;
  return new_position;
}

vector3d verlet_update(vector3d net_acceleration, vector3d previous_position)
{
  vector3d new_position = 2 * position - previous_position + delta_t * delta_t * net_acceleration;
  velocity = (new_position - position) / delta_t;
  return new_position;
}

void numerical_dynamics(void)
// This is the function that performs the numerical integration to update the
// lander's pose. The time step is delta_t (global variable).
{
  static vector3d previous_position;
  double fuel_mass = fuel * FUEL_CAPACITY * FUEL_DENSITY;
  double combined_mass = UNLOADED_LANDER_MASS + fuel_mass;

  vector3d gravity = -GRAVITY * MARS_MASS * combined_mass / position.abs2() * position.norm();
  vector3d drag = (parachute_status == DEPLOYED) ? drag_with_parachute(position) : drag_without_parachute(position);

  vector3d net_force = gravity + drag + thrust_wrt_world();
  vector3d net_acceleration = net_force / combined_mass;

  vector3d new_position = (simulation_time == 0.0) ? euler_update(net_acceleration) : verlet_update(net_acceleration, previous_position);

  previous_position = position;
  position = new_position;

  // Here we can apply an autopilot to adjust the thrust, parachute and attitude
  if (autopilot_enabled)
    autopilot();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude)
    attitude_stabilization(position.norm());
}

void initialize_simulation(void)
// Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "areostationry orbit";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";

  switch (scenario)
  {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2 * MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2 * MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE / 2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 6:
    // a areostationary equatorial orbit
    position = vector3d(AREOSTATIONARY_ORBIT_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, AREOSTATIONARY_ORBIT_SPEED, 0.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;
  }
}
