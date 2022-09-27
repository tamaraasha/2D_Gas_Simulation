# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 00:39:16 2019

@author: Tamara
"""
"""
simmod module
Tamara Anand 2019
"""

#%%
import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
import pylab as pl
from ballmod import *

class Simulation:
    
    def __init__(self, n): 
        """
        n initialises the simulation with n particles within a container
        At the default radii of container and particles, no more than 260 
        particles should be added to prevent overlap, n\<260.
        self._balls[0] should be the container. 
        """
        self._balls = [] #make first container, self._balls creates a list of each ball
        self._n = n
        self._time_elapsed = 0
        self._ketot = 0 #total KE of system
        self._meanv = 0
        self._ballmass= 1
        self._temp = 100
        self._pressuredeltav = [] #keeps track of impulse of each cont collision
        self._pressuretimes = [] #time of each contribution
        self._xmomovertime = []
        self._ymomovertime = []
        self._currentpr = [] #fill with improving mean pressures from 1 
                             #collision with container

            
    def kinetic(self):
        """
        kinetic() sums the kinetic energy values of each ball in self._balls
        and updates the value of total kinetic energy of the system, self._ketot
        """
        kes = [] 
        for i in range(0,len(self._balls)): 
            kes.append(self._balls[i]._ke) #for each of balls, add its KE to list
        sumkes = np.sum(kes) #sum all KEs
        self._ketot = sumkes #update system's KE
        return self._ketot
    
    

    def generate_balls(self, temp):
        """
        generate_balls fills simulation with a container and balls. The container's
        properties are hidden and set at default values. 
        Particles are given random velocities with a mean magnitude of the 
        velocity corresponding to the average velocity of particles in a system
        at equilibrium of the given temp.
        The velocities are randomly drawn from a uniform distribution whose 
        centre is this mean velocity.
        
        """
        self._ballrad = 0.1
        self._temp = temp
        self._contrad = 40 #container radius high so balls don't overlap
        self._contm = 1e10 #container mass high to minimise its movement
        self._contvel = np.array([0.0,0.0]) 
        self._contpos = np.array([0.0,0.0])
        ball1 = Ball(mass=self._contm, radius=self._contrad, position=self._contpos, vel=self._contvel, form=100)
       
        self._balls.append(ball1)
        pos = np.linspace(-39, 39, self._n) 
        kb = 1
      
        self._ballmass = 1 #all particles identical
        meanv = np.sqrt(2*kb*temp/self._ballmass) #1 = mass of ball
        self._meanv = meanv
        velx = rnd.uniform(0.0,2*meanv,self._n)
        vely = rnd.uniform(0.0,2*meanv,self._n)
        for i in range(pos.size):
            velx1 = velx[i]*rnd.uniform(-1,1)
            vely1 = vely[i]*rnd.uniform(-1,1)
            tempball = Ball(mass=self._ballmass, radius=self._ballrad, position=np.array([pos[i], 0.0]), vel=np.array([velx1, vely1]))
            self._balls.append(tempball)
        return self._balls
    

    def test_balls(self,ball):
        """
        Use Ball class to manually produce balls and add them to a simulation
        individually. 
        """
        self._balls.append(ball)
        self._n += 1
        return self._balls
    
    def set_part_rad(self,newrad):
        """
        Set the value of the particles' radii (excluding the container)
        """
        for i in range(1,len(self._balls)):
            self._balls[i]._rad = newrad
            

    def next_collision(self): 
        """
        Finds the next collision that will take place between any pair in the
        class instance and collides these particles, updating or maintaining
        the velocities and positions of all particles accordingly.
        """
        times = []
        collision_list = []
       
        for i, iball in enumerate(self._balls): #returns the index and 
            for j, jball in enumerate(self._balls): #object in self._balls array
                if j>i: #ensures not colliding with same ball 
                    t = self._balls[i].time_to_collision(self._balls[j]) 
                    """collision time for the ith ball with the jth ball"""
                    times.append(t) #add collision times to a list
                    collision_list.append((self._balls[i],self._balls[j])) 
                    """create a list of tuples containing the 'ID' of the 
                       colliding balls"""
        dt, idt = min(times), times.index(min(times)) 
        
        for i, balls in enumerate(self._balls): 
                 balls.move(0.99999*dt) 
                 #almost touching, if they actually touch, code doesn't progress
                                       
 
        ball1, ball2= collision_list[idt]  #identify colliding pair
                                          
      
        sbf1 = ball1._vel #collecting velocities before and after
        sbf2 = ball2._vel
       
        ball1.collide(ball2) 
       
        saf1 = ball1._vel
        saf2 = ball2._vel
        #note -  mass of particles is 1
        if ball2 == self._balls[0]: #if 1 collision ball is container
            deltavx = sbf1[0] - saf1[0] #changes in x and y momenta of ball
            deltavy = sbf1[1] - saf1[1]
            deltavtot = np.array([deltavx,deltavy]) 
            speeddif = np.linalg.norm(deltavtot) #magnitude of impulse
            self._pressuredeltav.append(speeddif) 
            self._pressuretimes.append(self._time_elapsed) 
        elif ball1 == self._balls[0]: 
            deltavx = sbf2[0] - saf2[0]
            deltavy = sbf2[1] - sbf2[1]
            deltavtot = np.array([deltavx,deltavy]) #if other ball is container
            speeddif = np.linalg.norm(deltavtot)
            self._pressuredeltav.append(speeddif)
            self._pressuretimes.append(self._time_elapsed)
        
        self._time_elapsed += dt
        return dt
   
    def pressure(self): #mean pressure from a simulation
        """
        num_frames in run function must be greater than 500 for a pressure
        to be generated for a simulation
        """
        pres = self._n*np.sum(self._currentpr[500:])/float(len(self._currentpr[500:]))
        return pres
 
    def run(self, num_frames, animate=False): 
        """
        num_frames equates to number of collisions to run. 
        If animate = True, animation of simulation is displayed.
        Plots of momentum, kinetic energy and particle distance from container
        centre over time are produced.
        """
        distances = []
        safter = []
        kes = []
        time = []
        sbinsta=[]
        vxinsta = []
        vyinsta = []
        allvy=[]
        mom_tot =[]

        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-40, 40), ylim=(-40, 40))
            ax.add_artist(self._balls[0].get_patch())
            for i in range (1,len(self._balls)):
                ax.add_patch(self._balls[i].get_patch())

    
        for frame in range(0,num_frames):
            vxinsta =[]
            vyinsta = []
            mom_frame=0.
            kes.append(self.kinetic()) #list of total KE over time
            time.append(self._time_elapsed) #list of corresponding times 
            pdv = self._pressuredeltav
            averagedeltav = np.sum(pdv)/len(pdv) #mean impulse from collisions
            circumf = 2*np.pi*40 #40 = container rad
            avpres = self._balls[1]._m*averagedeltav/(circumf*self._time_elapsed)
            self._currentpr.append(avpres)
            for i in range (1, len(self._balls)): #number of balls in container
                disti = np.linalg.norm(self._balls[i]._pos)
                vx = self._balls[i]._vel[0]
                vxinsta.append(vx)
                vy = self._balls[i]._vel[1]
                vyinsta.append(vy)
                allvy.append(vy)
                sb = np.linalg.norm(self._balls[i]._vel) #mag of each velocity
                sbinsta.append(sb) #adds each speed to array
                distances.append(disti)
                
            for i in range(0,len(self._balls)):
                mom_frame += np.dot(self._balls[i]._vel, 
                                    np.array([1.,0.]))*self._balls[i]._m
            mom_tot.append(mom_frame)
                
            self.next_collision()
            
            for i in range (1, len(self._balls)):
                 sa = np.linalg.norm(self._balls[i]._vel)
                 safter.append(sa)
        
            if animate:
                pl.pause(0.001)

            sbinsta=[] #sbinsta emptied after each frame
        if animate:
            pl.show()
            pl.xlim(-70,70)
        

        #Plot of Pressure
        plt.figure(1)
        plt.xlabel("Time (arbitrary units)")
        plt.ylabel("Pressure")
        plt.title("Pressure Over Time")
        plt.plot(self._currentpr, 'o') #Pressure with Time
         
        #Plot of x-component of Momentum
        plt.figure(1)
        xval = np.linspace(0,self._time_elapsed,2*num_frames)
        plt.xlabel('Time')
        plt.ylabel('Momentum')
        plt.title('Total Momentum Over Time')
        plt.plot(mom_tot)
        
        #Histogram of Particle Positions
        distarr = np.asarray(distances)  
        plt.figure(1)
        plt.xlabel('Distance from centre')
        plt.ylabel('Frequency')
        plt.title('Histogram of Distance of Balls from Container Centre (R = 40)')
        plt.hist(distarr)
        
        mean = np.mean(safter)
        meanan = "Mean Speed is", mean
        def MB_dis(v):
            multby = v*self._ballmass/(1*self._temp) #where kB is taken as 1
            exp = np.exp(-(v**2)*self._ballmass/(2*1*self._temp))
            y = multby*exp
            return y
        ss = np.arange(0,100,0.1)
        
        meanvy = np.mean(allvy)
        meananvy = "Mean Velocity is", meanvy
        print(self._temp)
        def MB_disv(v):
            multby = np.sqrt(self._ballmass/(2*np.pi*self._temp)) #where kB is 
            exp = np.exp(-(v**2)*self._ballmass/(2*1*self._temp)) #taken as 1
            y = multby*exp
            return y
        
        
        def MB_var(temp): #calculates variance of MB velocity distribution
            asqr = temp #would equal kb*T/m but other values = 1
            var = asqr*(3*np.pi - 8)/np.pi
            return var
        
        #Histogram of y-component velcities with MB distribution
        plt.figure(4) 
        vs = np.arange(-100,100,0.1)
        plt.plot(vs,MB_disv(vs))
        plt.xlabel('Velocity of Ball') 
        plt.ylabel('Frequency')
        plt.title("Histogram showing Velocity Distribution (Vertical Component)") 
        plt.hist(allvy,bins=85,density=True)
        
        #Histogram of speeds with MB distribution        
        plt.figure(5) 
        plt.plot(ss,MB_dis(ss))
        plt.xlabel('Speed of Ball') 
        plt.ylabel('Frequency')
        plt.title('Histogram showing Ball Speed Distribution') 
        plt.hist(safter,bins=85,density=True)
        variancesim = np.var(safter) 
        variancespred = MB_var(self._temp)
        print("Variance - Simulated Data", variancesim, 
              "Versus Variance - MB fit", variancespred)
        
        #Plot of kinetic energy over time
        plt.figure(6) 
        plt.xlabel('Time Elapsed')  
        plt.ylabel('Kinetic Energy of the System')
        plt.title('Plot of Kinetic Energy of the System over Time')
        plt.plot(time,kes,'o')
