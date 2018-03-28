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
    global jup, sun # Read in the celestial objects
    global dist # Read in the solar-Jupiter distance
    r = np.linalg.norm(jup.r) # Get the radius of Jupiter-like
    jup_angle = np.arcsin(jup.r0[1]/r) # Initial angle of Jupiter-like
    planet_angle = jup_angle + angle
    planet_x = dist*np.cos(planet_angle)+sun.r0[0] # Calculate the angle from Sun position!
    planet_y = dist*np.sin(planet_angle)+sun.r0[1]
    planet_pos = [planet_x+dr[0], planet_y+dr[1], dr[2]] # Unperturbed pos is in plane z=0
    # Then also calculate the "circular velocity" of steady orbit at same speed as Jup
    planet_vx = -2*np.pi*np.sqrt(r/dist**2)*np.sin(planet_angle)
    planet_vy = 2*np.pi*np.sqrt(r/dist**2)*np.cos(planet_angle)
    planet_vel = [planet_vx+dv[0], planet_vy+dv[1], dv[2]] # Add velocity perturbation
    # Now perturbed position and circular velocities calculated, return 6-vector:
    return planet_pos+planet_vel

def plot_solution_xy(t_points, soln, fig_num):
    # Function to plot both Trojan and Jupiter-like solutions on the same x-y plane map.
    # First get the Jupiter/Su solutions at the time points:
    sunjup_soln = exact_soln(t_points)
    sun_pos = np.array(exact_soln(t_points)[1]) # Make both solutions into arrays
    jup_pos = np.array(exact_soln(t_points)[0]) 
    troj_pos = np.array(soln) # To make for easy indexing
    plt.figure(fig_num) # Make the figure
    plt.plot(jup_pos[:,0],jup_pos[:,1],'b-')
    # plt.plot(sun_pos[:,0],sun_pos[:,1],'gx')
    plt.plot(troj_pos[:,0],troj_pos[:,1],'r-')
    plt.show()

def trojan_stability_measure(t_points, soln, sign):
    # Function which takes the end point of the asteroid at the final time point
    # provided, and then calculates the exact Jupiter point at that time, and then
    # the theoretical exact solution for the final asteroid position (given that
    # it was at sgn(sign)*pi/3 from Jupiter-like), then returns the difference as
    # a fraction of the Jupiter radius. 
    global dist # Read in the Sun-Jupiter distance
    final_t = [t_points[-1]] # Select final slice for exact solution
    jup_soln = np.array(exact_soln(final_t)[0]+[0]) # Gives the Jup position at final t_point
    sun_soln = np.array(exact_soln(final_t)[1]+[0]) # Gives the Sun position at final t_point
    trojan_soln = np.array(soln[-1][0:3]) # Gives the calculated Trojan position at final t_points
    jup_r = np.linalg.norm(jup_soln)
    jup_angle = np.arctan2(jup_soln[1],jup_soln[0]) # Unambiguous angle from -pi > pi
    trojan_angle = jup_angle + sign*(np.pi/3) # Generate angle planet would be at
    theor_trojan_soln = sun_soln + np.array([dist*np.cos(trojan_angle),dist*np.sin(trojan_angle),0])
    trojan_diff_dist = np.linalg.norm(trojan_soln-theor_trojan_soln)
    trojan_frac_diff = trojan_diff_dist/jup_r
    return trojan_frac_diff


angle = np.pi/3
dr = (0*np.array([np.cos(angle),np.sin(angle),0])).tolist()
dv = (0.1*np.array([0,0,1])).tolist()
y0 = perturbed_initial_condition(dr,dv,angle)
trojan = Planet(y0,0) # Create the Trojan object (mass unimportant)
t_points = np.linspace(0,1000*11.8,100000)
soln = integrate_gravity(y0, t_points, [jup,sun]) 
plot_solution_xy(t_points, soln, 1)
# print(trojan_stability_measure(t_points, soln, -1))
# plt.figure(2)
# asoln = np.array(soln)
# plt.plot(asoln[:,2],'b-')
# plt.show()
