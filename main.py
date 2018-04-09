# File which uses an adaptive LSODA method in order to investigate the orbits of
# negligible masses around the Lagrange points of large mass planet-star systems, to 
# work out the scaling relationship of the stable orbit region with the planetary mass.
# Units used in the problem are length unit 1AU, time unit 1yr, mass unit 1solarmass
# In this file, the frame of reference is the inertial frame of the Solar System,
# so the Sun, Jupiter and Lagrange point all rotate in this frame, no inertial forces.

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
jup_m = 0.001 # Control the mass of the Jupiter-like for the scaling relation
jup = Planet([0.0,0.0,0.0,0.0,0.0,0.0],jup_m) # Create the Jupiter-like (positions not defined)
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
    init_angle_jup = np.arctan2(jup.r0[1],jup.r0[0]) # Using unambiguous angle
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
    init_angle_sun = init_angle_jup + np.pi # Always starts pi away from Jupiter posn
    sun_r = np.linalg.norm(sun.r) # Get rotational radius of the Sun
    solar_positions = []
    for t_point in t:
        # Go through each time point
        # Then use the angular velocity to add on the angle for a given time t/yr since start
        new_angle = init_angle_sun + (angular*t_point)
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

def trojan_stability_measure(t_points, soln, sign):
    # Function which takes the end point of the asteroid at the final time point
    # provided, and then calculates the exact Jupiter point at that time, and then
    # the theoretical exact position of the Lagrance point being used (given that
    # it was at sgn(sign)*pi/3 from Jupiter-like), then returns the difference as
    # a fraction of the Jupiter radius. 
    global dist # Read in the Sun-Jupiter distance
    final_t = [t_points[-1]] # Select final slice for exact solution
    total_exact_soln = exact_soln(final_t)
    jup_soln = np.array(total_exact_soln[0]+[0]) # Gives the Jup position at final t_point
    sun_soln = np.array(total_exact_soln[1]+[0]) # Gives the Sun position at final t_point
    trojan_soln = np.array(soln[-1][0:3]) # Gives the calculated Trojan position at final t_points
    jup_r = np.linalg.norm(jup_soln)
    jup_angle = np.arctan2(jup_soln[1],jup_soln[0]) # Unambiguous angle from -pi > pi
    trojan_angle = jup_angle + sign*(np.pi/3) # Generate angle planet would be at
    theor_trojan_soln = sun_soln + np.array([dist*np.cos(trojan_angle),dist*np.sin(trojan_angle),0])
    trojan_diff_dist = np.linalg.norm(trojan_soln-theor_trojan_soln)
    trojan_frac_diff = trojan_diff_dist/jup_r
    return trojan_frac_diff

def perturbed_initial_condition(dr,dv,angle):
    # Function which produces the initial condition 6-vector for a planet
    # y0 = [x,y,z,vx,vy,vz] based on its radian  angle from Jupiter (negative is counted
    # as clockwise from Jupiter) and a perturbation (dx,dy,dz) from that point.
    # The initial velocity is given by the 3-vector v.
    global jup, sun # Read in the celestial objects
    global dist # Read in the solar-Jupiter distance
    jup_r = np.linalg.norm(jup.r) # Get the radius of Jupiter-like
    jup_angle = np.arcsin(jup.r0[1]/jup_r) # Initial angle of Jupiter-like
    new_angle = jup_angle + angle # Angle wrt the Sun
    planet_x = dist*np.cos(new_angle)+sun.r0[0] # Calculate the angle from Sun position!
    planet_y = dist*np.sin(new_angle)+sun.r0[1]
    planet_pos = [planet_x+dr[0], planet_y+dr[1], dr[2]] # Unperturbed pos is in plane z=0
    planet_r = np.linalg.norm([planet_x,planet_y])
    # Then also calculate the "circular velocity" of orbit to match Jupiter angular speed
    trojan_angle = np.arctan2(planet_y,planet_x) # This is the undisturbed angle wrt the ORIGIN
    planet_vx = -2*np.pi*(planet_r/dist)*np.sqrt(1/jup_r)*np.sin(trojan_angle)
    planet_vy = 2*np.pi*(planet_r/dist)*np.sqrt(1/jup_r)*np.cos(trojan_angle)
    planet_vel = [planet_vx+dv[0], planet_vy+dv[1], dv[2]] # Add velocity perturbation
    # Now perturbed position and circular velocities calculated, return 6-vector:
    return planet_pos+planet_vel
