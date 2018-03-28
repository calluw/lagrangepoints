# File which uses an adaptive Runge-Kutta method in order to investigate the orbits of
# negligible masses around the Lagrange points of large mass planet-star systems, to 
# work out the scaling relationship of the stable orbit region with the planetary mass.
# Units used in the problem are length unit 1AU, time unit 1yr, mass unit 1solarmass

## Imports
import numpy as np
from scipy.integrate import odeint

class Planet:
    def __init__(self,y0,m):
        # Create the planet/celestial object with its position r and velocity v as 3D vectors
        # Mass m defined, for test masses set m=0, for Jupiter mass/Solar mass use m/=0
        self.r0 = np.array([y0[0],y0[1],y0[2]])
        self.r = self.r0
        self.v0 = np.array([y0[3],y0[4],y0[5]])
        self.v = self.v0
        self.m = m

## Globals
jup = Planet([0.0,0.0,0.0,0.0,0.0,0.0],0.001) # Create the Jupiter-like m=0.001 (positions not defined)
sun = Planet([0.0,0.0,0.0,0.0,0.0,0.0],1) # Create a Solar mass
dist = 5.2 # Define the average distance between the 2 bodies

def update_barycentre():
    # Define the barycentre of the Sun-Jupiter-like system as [0,0,0]
    # Hence using M_j/M_sol, and their distance d, get distance of Jup/Sun:
    global jup, sun # Read in Jupiter-like and Solar objects to update
    global dist # Read in their average separation
    rad_s = dist/(1+(sun.m/jup.m)) # Definition of barycentre
    rad_j = dist-rad_s
    sun.r0[0] = -rad_s # Create a Sun in -x dirn
    jup.r0[0] = rad_j # And a Jupiter in +x dirn
    sun.r = sun.r0
    jup.r = jup.r0 # Update their positions just in case

update_barycentre() # Ensure the barycentre distances are updated for future calculations

def exact_soln(t):
    # Function which uses the planet object for Jupiter-like jup and the solar object
    # sun in order to construct the known analytic solution for the positions given
    # a time from initial position, under the assumption that corotational orbit is in x-y
    # This function can accept either a single time point in [] or an array of them as t.
    
    global jup, sun # Import the celestial objects 
    global dist # Distance needed for angular velocity calc
    jup_r = np.linalg.norm(jup.r) # Get rotational radius of Jupiter-like
    angular = 2*np.pi*np.sqrt(1/(jup_r))*(1/dist) # A constant angular velocity in rad/yr for both
    
    ## Jupiter
    # Next calculate initial angle, angle defined as 0 radians being Jupiter at x=r,y=0
    init_angle_jup = np.arcsin(jup.r0[1]/jup_r) # Using sin(theta) = opp/hyp
    jupiter_positions = []
    for t_point in t:
        # Go through each time point
        # Then use the angular velocity to add on the angle for a given time t/yr since start
        new_angle = init_angle_jup + (angular*t_point)
        jup_x = jup_r*np.cos(new_angle)
        jup_y = jup_r*np.sin(new_angle)
        jupiter_positions.append([jup_x,jup_y]) 
    if len(jupiter_positions) < 2:
        jupiter_positions = jupiter_positions[0] # If length 1, remove the outer list layer 

    ## Sun
    init_angle = init_angle_jup + np.pi # Always starts pi away from Jupiter posn
    sun_r = np.linalg.norm(sun.r) # Get rotational radius of the Sun
    solar_positions = []
    for t_point in t:
        # Go through each time point
        # Then use the angular velocity to add on the angle for a given time t/yr since start
        new_angle = init_angle + (angular*t_point)
        sun_x = sun_r*np.cos(new_angle)
        sun_y = sun_r*np.sin(new_angle)
        solar_positions.append([sun_x,sun_y]) 
    if len(solar_positions) < 2:
        solar_positions = solar_positions[0] # If length 1, remove the outer list layer 

    return [jupiter_positions, solar_positions]


def grav_derivatives(y, t, objects):
    # Function which calculates the derivatives of the planet's position and velocity
    # at a given time t, using the gravitational forces from the objects provided.
    # Note that the 6 coordinates are defined with 3 position in self.r, 3 velocities in self.v.
    
    r = np.array([y[0],y[1],y[2]]) # Create the self position vector
    v = np.array([y[3],y[4],y[5]]) # Create the self velocity vector
    r_dash = v # Trivial ODEs from time derivates
    exact = exact_soln([t]) # Get the exact solution at t
    jup.r = np.array(exact[0]+[0]) # Update Jupiter position in plane z=0 at single time point t
    sun.r = np.array(exact[1]+[0]) # Update Solar position in plane z=0 at single time point t
    grav_acc = np.zeros(3) # Make acc vector, note that acceleration will be in (AU s^-2)
    for object in objects:
        dist_vec = r-object.r # Create the vectorial distance
        dist_to_obj = np.linalg.norm(dist_vec) # And the scalar distance
        f_strength = (-4*(np.pi**2)*object.m)/(dist_to_obj**3) # The 'strength' of the force
        grav_acc += f_strength*dist_vec # Acceleration += (-GM/r^3)*_r_, update acceleration
    v_dash = grav_acc # Derivative of velocities are the accelerations
    return r_dash.tolist()+v_dash.tolist() # Return the derivatives as concatenated 6-vector

def integrate_gravity(y0, t_points, objects):
    # Function which takes the initial conditions y0 as 6-vector [r1,r2,r3,v1,v2,v3] for 3D space
    # and a series of time points at which to evaluate the solutions it finds (via adaptive
    # stepping) given a landscape of massive bodies in the objects array of Planet class objects.
    # Outputs the y 6-vector at the time points t_points.
    
    y = odeint(grav_derivatives,y0,t_points,args=(objects,)) # Calling the integrator function
    return y # Return a list giving the y 6-vector at each time point

