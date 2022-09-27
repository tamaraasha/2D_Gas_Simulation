# thermo_simulation

Monatomic Gas Simulation - Tamara Anand 2019

This simulation consists of 2 modules each corresponding to a class.

1.	ballmod

The module ballmod contains the class "Ball", which produces individual balls of a given mass, radius, position and velocity. 
It also contains functions involving the collision of 2 balls. time_to_collision finds the time until a collision will take
place between two balls, and collide updates the velocities and kinetic energies of the balls, effectively performing an 
instantaneous collision. 

2. 	simmod

This module contains the class "Simulation". Simulation builds on the previous class in that it can use balls that Ball 
creates and store collections of them within a system. generate_balls takes temperature as an input to produce n balls
whose mean velocities correspond to what would be expected of a gas at equilibrium at the given temperature. test_balls
allows the addition of manually designed balls to any system. The run function finally runs a given simulation for a given number of frames, with each frame corresponding to a collision. It automatically produces graphs showing the evolution of pressure, total momentum and total kinetic energy over time. Histograms of speed, velocity and position for each ball over every frame are also produced.

The scripts file imports both classes and contains a number of cells, each of which have examples that correspond to various tasks. 

List of files submitted:
- simmod.py
- ballmod.py
- scripts.py
