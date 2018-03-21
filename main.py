# File which uses an adaptive Runge-Kutta method in order to investigate the orbits of
# negligible masses around the Lagrange points of large mass planet-star systems, to 
# work out the scaling relationship of the stable orbit region with the planetary mass.
# Units used in the problem are length unit 1AU, time unit 1yr, mass unit 1solarmass

import numpy as np
import matplotlib.pyplot as plt

class Planet:
    def __init__(self,rx,ry,rz,vx,vy,vz,m):
        # Define the planet/celestial object with its position r and velocity v as 3D vectors
        # Mass m defined, for test masses set m=0, for Jupiter mass/Solar mass use m/=0
        self.r0 = np.array([rx,ry,rz])
        self.r = self.r0
        self.v0 = np.array([vx,vy,vz])
        self.v = self.v0
        self.m = m

def jup_soln(jup,t):
    # Function which uses the planet object for Jupiter-like jup and the solar object
    # sun in order to construct the known analytic solution for the position of jup given
    # a time from its initial position, under the assumption that jup orbit is in x-y plane
    # The analytic solution in these units gives v_jup = 2pi*sqrt(1/r), and R_jup is 5.2AU 
    # (5.2 units), hence the angular velocity v=rw => w = 2pi*sqrt(1/r^3)
    angular = 2*np.pi*np.sqrt(1/5.2**3) # A constant in rad/yr
    # Next calculate initial angle, angle defined as 0 radians being Jupiter at x=5.2,y=0
    init_angle = np.arcsin(jup.r0[1]/5.2) # Using sin(theta) = opp/hyp
    print(init_angle)
    # Then use the angular velocity to add on the angle for a given time t/yr since start
    new_angle = init_angle + (angular*t)
    print(new_angle)
    jup_x = 5.2*np.cos(new_angle)
    jup_y = 5.2*np.sin(new_angle)
    return [jup_x,jup_y]

jup = Planet(5.2,0,0,0,0,0,0.001)
print(jup_soln(jup,1))
# Program currently doesn't work as there is an issue with the exact solution (v or w likely wrong)
