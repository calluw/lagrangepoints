# File which takes the code from main.py and the uses it to investigate the orbits 
# described there. Essentially this is an implementation file which records how 
# the output of the integration in main.py was used to produce images and observations
# used in the report.

## Imports
import numpy as np
import matplotlib.pyplot as plt
from main import * # Import all functionality and globals from the main.py file

def perturbedInitialCondition(dr,v,angle):
    # Function which produces the initial condition 6-vector for a planet
    # y0 = [x,y,z,vx,vy,vz] based on its radian  angle from Jupiter (negative is counted
    # as clockwise from Jupiter) and a perturbation (dx,dy,dz) from that point.
    # The initial velocity is given by the 3-vector v.
    global jup # Read in the Jupiter-like for its initial position
    jup_angle = np.arcsin(jup.r0[1]/5.2) # Initial angle of Jupiter-like
    planet_angle = jup_angle + angle
    planet_x = 5.2*np.cos(planet_angle)
    planet_y = 5.2*np.sin(planet_angle)
    planet_pos = [planet_x+dr[0], planet_y+dr[1], dr[2]] # Unperturbed pos is in plane z=0
    # Now perturbed position calculated, return 6-vector:
    return planet_pos+v

y0 = perturbedInitialCondition([0,0,0],[0,0,0],-np.pi/3) # Initial unperturbed Trojan at rest
trojan = Planet(y0,0) # Create the Trojan object (mass unimportant)
t_points = np.linspace(0,2,1000)
y = integrate_gravity(y0, t_points, [jup,sun])
for line in y:
    print(line)
    
