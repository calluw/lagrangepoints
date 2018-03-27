# File which takes the code from main.py and the uses it to investigate the orbits 
# described there. Essentially this is an implementation file which records how 
# the output of the integration in main.py was used to produce images and observations
# used in the report.

## Imports
import numpy as np
import matplotlib.pyplot as plt
from main import * # Import all functionality and globals from the main.py file

def perturbed_initial_condition(dr,dv,angle):
    # Function which produces the initial condition 6-vector for a planet
    # y0 = [x,y,z,vx,vy,vz] based on its radian  angle from Jupiter (negative is counted
    # as clockwise from Jupiter) and a perturbation (dx,dy,dz) from that point.
    # The initial velocity is given by the 3-vector v.
    global jup # Read in the Jupiter-like for its initial position
    r = np.linalg.norm(jup.r) # Get the radius of Jupiter-like
    jup_angle = np.arcsin(jup.r0[1]/r) # Initial angle of Jupiter-like
    planet_angle = jup_angle + angle
    planet_x = r*np.cos(planet_angle)
    planet_y = r*np.sin(planet_angle)
    planet_pos = [planet_x+dr[0], planet_y+dr[1], dr[2]] # Unperturbed pos is in plane z=0
    # Then also calculate the "circular velocity" of steady orbit at same speed as Jup
    planet_vx = -2*np.pi*np.sqrt(1/r)*np.sin(planet_angle)
    planet_vy = 2*np.pi*np.sqrt(1/r)*np.cos(planet_angle)
    planet_vel = [planet_vx+dv[0], planet_vy+dv[1], dv[2]] # Add velocity perturbation
    # Now perturbed position and circular velocities calculated, return 6-vector:
    return planet_pos+planet_vel

def plot_solution_jup(t_points, soln, fig_num):
    # Function to plot both Trojan and Jupiter-like solutions on the same x-y plane map.
    # First get the Jupiter solution at the time points:
    jup_pos = np.array(jup_soln(t_points)) # Make both solutions into arrays
    troj_pos = np.array(soln) # To make for easy indexing
    plt.figure(fig_num) # Make the figure
    #  plt.plot(jup_pos[:,0],jup_pos[:,1],'b-')
    plt.plot(troj_pos[:,0],troj_pos[:,1],'r-')
    plt.show()


# y0 = perturbed_initial_condition([0,0,0],[0,0,0],np.pi/3) # Initial Trojan 
# trojan = Planet(y0,0) # Create the Trojan object (mass unimportant) with circular orbital velocity
# t_points = np.linspace(0,11800,20000)
# soln = integrate_gravity(y0, t_points, [jup,sun]) 
# plot_solution_jup(t_points, soln, 1)

