# File which takes the code from main.py and the uses it to investigate the orbits 
# described there. Essentially this is an implementation file which records how 
# the output of the integration in main.py was used to produce images and observations
# used in the report.

## Imports
import numpy as np
import matplotlib.pyplot as plt
from main import * # Import all functionality and globals from the main.py file
from main_rotatingframe import * # Including the rotating frame variants

def plot_solution_xy(t_points, soln, fig_num):
    # Function to plot both Trojan and Jupiter-like solutions on the same x-y plane map.
    # First get the Jupiter/Su solutions at the time points:
    sunjup_soln = exact_soln(t_points)
    sun_pos = np.array(exact_soln(t_points)[1]) # Make both solutions into arrays
    jup_pos = np.array(exact_soln(t_points)[0]) 
    troj_pos = np.array(soln) # To make for easy indexing
    plt.figure(fig_num) # Make the figure
    # plt.plot(jup_pos[:,0],jup_pos[:,1],'b-')
    # plt.plot(sun_pos[:,0],sun_pos[:,1],'gx')
    plt.plot(troj_pos[:,0],troj_pos[:,1],'r-')
    # plt.show()

angle = -np.pi/3 # Leading or trailing group (+ve is leading)
dr = (0*np.array([np.cos(angle),np.sin(angle),0])).tolist() # Generate positional perturbation
dv = (0*np.array([-np.sin(angle),np.cos(angle),0])).tolist() # Generate velocity perturbation
# y0 = perturbed_initial_condition(dr,dv,angle) # Initial position in the case of inertial frame
y0 = perturbed_initial_condition_rot(dr,dv,angle) # Initial position in case of rotating frame
trojan = Planet(y0,0) # Create the Trojan object (mass unimportant)
t_points = np.linspace(0,500*11.8,20000) # Define the time points used in simulation
# soln = integrate_gravity(y0, t_points, [jup,sun]) # Perform the integration (inertially)
soln = integrate_gravity_rot(y0, t_points, [jup,sun]) # Integrate in rotating frame (inertial forces)

plot_solution_xy(t_points, soln, 1) # Plot solution in x-y plane
# print(trojan_stability_measure(t_points, soln, 1)) # Give the measure of the stability of soln
print(trojan_stability_measure_rot(t_points, soln, -1))
# plt.figure(2)
# total_fraction_arr = [] # Get the fractional distance to show lateral displacements
# for i,t in enumerate(t_points):
#     total_fraction_arr += [trojan_stability_measure_rot([t], soln[i:i+1], +1)]
# print(total_fraction_arr[0]) # Print the initial displacement fraction
# print(total_fraction_arr[-1])
# plt.plot(total_fraction_arr)
plt.show()

