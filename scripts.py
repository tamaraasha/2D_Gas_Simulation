# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 00:40:52 2019

@author: Tamara
"""

import numpy as np
import matplotlib.pyplot as plt
import ballmod as ba
import simmod as sm
#%% TASK 4 add balls manually
s = sm.Simulation(0)
b1 = ba.Ball(mass=1, radius=1, position=np.array([-5.0, 0.0]), 
             vel=np.array([1.0, 0.0]), form=10)
b0 = ba.Ball(mass=100000000, radius=20, position=np.array([0.0, 0.0]),
             vel=np.array([0.0, 0.0]), form=100)
s.test_balls(b0)
s.test_balls(b1)
s.run(100,animate = True)
#%% VEOCITY IN TWO DIRECTIONS
s = sm.Simulation(0)
b3 = ba.Ball(mass=1, radius=1, position=np.array([-5.0, 7.0]), 
             vel=np.array([1.0, 0.0]), form=10)
b2 = ba.Ball(mass=100000000, radius=20, position=np.array([0.0, 0.0]), 
             vel=np.array([0.0, 0.0]), form=100)
s.test_balls(b2)
s.test_balls(b3)

s.run(100,animate = True)

#%% EXAMPLE SIMULATION AT 400K
exm = sm.Simulation(30)
exm.generate_balls(400) 
exm.run(30000,animate=False)
print("Average Pressure",exm.pressure())
print("Total Kinetic Energy",exm.kinetic())
#%% SHOWING DIFFERENT MAXWELL BOLTZMANN DISTRIBUTIONS AT DIFFERENT TEMPERATURES
x = sm.Simulation(35)
x.generate_balls(100)
y = sm.Simulation(35)
y.generate_balls(500)
z = sm.Simulation(35)
z.generate_balls(900)
x.run(30000, animate = False)
y.run(30000, animate = False)
z.run(30000, animate = False)
#%% TASK 12, RELATION BETWEEN PRESSURE AND TEMPERATURE
on = sm.Simulation(30) #had attempted production of these over for loop but
on.generate_balls(100) #was unsuccesful, 6 data points is enough to indicate
tw = sm.Simulation(30) #trend in any case
tw.generate_balls(500)
th = sm.Simulation(30)
th.generate_balls(900)
fo = sm.Simulation(30)
fo.generate_balls(1300)
fi = sm.Simulation(30)
fi.generate_balls(1700)
si = sm.Simulation(30)
si.generate_balls(2100)
on.run(10000, animate = False)
tw.run(10000, animate = False)
th.run(10000, animate = False)  
fo.run(10000, animate = False)
fi.run(10000, animate = False)
si.run(10000, animate = False)
pressure = [on.pressure(),tw.pressure(),th.pressure(),fo.pressure(),
            fi.pressure(),si.pressure()]
t = [100,500,900,1300,1700,2100]

idealgrad = 1/(np.pi*(40)**2) #this is the predicted gradient from ideal gas eqn, 40 is container radius

temps = np.linspace(0,2500,1000)
a,b = np.polyfit(t, pressure,1) 
def fitt(x,a,b): #best fit equation for simulated data
    y = a*x + b
    return y

plt.figure(10)
plt.plot(temps,idealgrad*temps) #predicted plot for ideal gas
plt.plot(temps,fitt(temps,a,b))
plt.plot(t,pressure,'o')
plt.xlabel("Absolute Temperature (K)")
plt.ylabel("Pressure")
plt.title("Plot of Pressure vs Temperature")
print("IDEAL GAS GRADIENT", idealgrad, "VS SIMULATION GRADIENT",a)
# =============================================================================
# 
# Doesn't seem to fit trend expected when particle radius is changed manually
# with self._ballrad for each simulation. 
# 
# =============================================================================
#%% START OF BROWNIAN MOTION INVESTIGATION
exm = sm.Simulation(80)
bigball = ba.Ball(mass=15, radius=1, position=np.array([0.0, 2.0]), vel=np.array([0.0001, 0.0]), form=10)
big = exm.test_balls(bigball)
exm.generate_balls(400) 
exm.run(2000,animate=True)
"""
I didn't have time to develop this investigation, but by adding a larger, 
heavier ball (out of the axis of particles), Brownian motion could be
investigated. A higher number of particles is important for the animation
to demonstrate this random motion of the larger ball. 
"""