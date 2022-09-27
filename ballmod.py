# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 00:38:37 2019

@author: Tamara
"""
"""
ballmod module
Tamara Anand 2019
"""

#%%
import numpy as np
import numpy.random as rnd
import pylab as pl
import matplotlib.pyplot as plt

class Ball:
    def __init__(self, mass=1e20, radius=10, position=np.array([0.0, 0.0]),
                 vel=np.array([0.0, 0.0]), form=1): 
        """
        Set form = 100 if the object is a container
        Mass of container must be many orders of magnitude bigger than particle
        mass
        """
        self._pos = position
        self._vel = vel 
        self._rad = radius
        self._m = mass 
        self._form = form 
        self._ke = 0.5*self._m*(np.linalg.norm(self._vel))**2 
     
        if form == 100:
            self._patch = pl.Circle(self._pos, self._rad, fc='b', fill=False) 
        else:
            self._patch = pl.Circle(self._pos, self._rad, fc='r')
       
    def __repr__(self):
        """
        Returns key properties of ball instance
        """
        return "%s,%g,%g,%s,%s" % ("Mass, Radius, Vel, Pos:", self._m,
                                   self._rad, self._vel, self._pos)
    
    def pos(self):
        return self._pos
    
    def vel(self):
        return self._vel
    
    def get_patch(self):
        return self._patch 
       
    def move(self, dt): #set function
        self._pos += dt*self._vel
        return self._pos
       
    def time_to_collision(self, other):
        """
        time_to_collision calculates solutions to 2 versions of a quadratic equation,
        producing 4 possible solutions to the time to next collision between
        2 particles. The returned value depends on if a container is involved 
        and the exact values produced.
        If no solution to quadratic and balls never collide, dt = 1e10 
        so collisions that will actually take place are lower and returned instead.
        Otherwise, the minimum positive value is returned
        """
        r_rel = self._pos - other._pos
        v_rel = self._vel - other._vel

        rad_add = self._rad + other._rad #if 2 particles colliding
        rad_dif = self._rad - other._rad #if colliding with container
        dot = np.dot(r_rel,v_rel) 
        vsquared = np.dot(v_rel,v_rel)
        rsquared = np.dot(r_rel,r_rel) 
        if self._form == 100 or other._form == 100: #run through container quadratic equation
            time1 = (-dot + np.sqrt(dot**2 - vsquared*(rsquared-rad_dif**2)))/vsquared
            time2 = (-dot - np.sqrt(dot**2 - vsquared*(rsquared-rad_dif**2)))/vsquared
            if time1<0.0 and time2>0.0: 
                dt = time2
                return time2
            elif time2<0.0 and time1>0.0:
                dt = time1
                return time1
            elif time1<0.0 and time2<0.0: #no collision
                dt = 1e10
                return 1e10 
            else:
                return min(time1,time2)
                dt = min(time1,time2)
        else: #2 balls colliding

            time3 = (-dot + np.sqrt(dot*dot - vsquared*(rsquared-rad_add**2)))/vsquared
            time4 = (-dot - np.sqrt(dot*dot - vsquared*(rsquared-rad_add**2)))/vsquared

            if time3<0.0 and time4>0.0:
                dt = time4
                return time4
            elif time4<0.0 and time3>0.0:
                dt = time3
                return time3
                
            elif time3<0.0 and time4<0.0: #no collision
                dt = 1e10
                return 1e10 

            else:
                return min(time3,time4)
                dt = min(time3,time4)
        
            
            
    def collide(self, other):
        """
        collide function considers the masses and velocities of two objects
        that will collide and calculates the new velocities based on the 2D
        equations for momentum conservation in elastic collisions.
        Velocities and kinetic energies are updated after collision. 
        """
        massadd = self._m + other._m
        r_rel = self._pos - other._pos
        v_rel = self._vel - other._vel
        magr = np.linalg.norm(r_rel)
        self._vel= self._vel - (2*other._m*np.dot(v_rel,r_rel)*r_rel)/(massadd*(magr)**2)
        other._vel= other._vel - (2*self._m*np.dot(-v_rel,-r_rel)*(-r_rel))/(massadd*(magr)**2) #this doesnt work if containter and ball centered on same point
        self._ke = 0.5*self._m*np.linalg.norm(self._vel)**2
        other._ke = 0.5*other._m*np.linalg.norm(other._vel)**2

        return self._vel, other._vel

