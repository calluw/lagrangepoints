# File which uses an adaptive Runge-Kutta method in order to investigate the orbits of
# negligible masses around the Lagrange points of large mass planet-star systems, to 
# work out the scaling relationship of the stable orbit region with the planetary mass.

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

 def jupiter_solution(jup,sun):
    # Function which uses the planet object for Jupiter-like jup and the solar object
    # sun in order to construct the known analytic solution for the position of jup given
    # a time from its initial position

