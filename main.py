# File which uses an adaptive Runge-Kutta method in order to investigate the orbits of
# negligible masses around the Lagrange points of large mass planet-star systems, to 
# work out the scaling relationship of the stable orbit region with the planetary mass.
# Units used in the problem are length unit 1AU, time unit 1yr, mass unit 1solarmass

## Imports
import numpy as np
from scipy.integrate import odeint

class Planet:
    def __init__(self,rx,ry,rz,vx,vy,vz,m):
        # Create the planet/celestial object with its position r and velocity v as 3D vectors
        # Mass m defined, for test masses set m=0, for Jupiter mass/Solar mass use m/=0
        self.r0 = np.array([rx,ry,rz])
        self.r = self.r0
        self.v0 = np.array([vx,vy,vz])
        self.v = self.v0
        self.m = m

## Globals
jup = Planet(5.2,0,0,0,0,0,0.001) # Create the Jupiter-like m=0.001 at x=5.2 (v is unimportant)
sun = Planet(0,0,0,0,0,0,1) # Create a Solar mass at the origin, which does not move

def jup_soln(t):
    # Function which uses the planet object for Jupiter-like jup and the solar object
    # sun in order to construct the known analytic solution for the position of jup given
    # a time from its initial position, under the assumption that jup orbit is in x-y plane
    # The analytic solution in these units gives v_jup = 2pi*sqrt(1/r), and R_jup is 5.2AU 
    # (5.2 units), hence the angular velocity v=rw => w = 2pi*sqrt(1/r^3)
    
    global jup # Import the jupiter object 
    angular = 2*np.pi*np.sqrt(1/5.2**3) # A constant in rad/yr
    # Next calculate initial angle, angle defined as 0 radians being Jupiter at x=5.2,y=0
    init_angle = np.arcsin(jup.r0[1]/5.2) # Using sin(theta) = opp/hyp
    # Then use the angular velocity to add on the angle for a given time t/yr since start
    new_angle = init_angle + (angular*t)
    jup_x = 5.2*np.cos(new_angle)
    jup_y = 5.2*np.sin(new_angle)
    return [jup_x,jup_y]


def grav_derivatives(y, t, objects):
    # Function which calculates the derivatives of the planet's position and velocity
    # at a given time t, using the gravitational forces from the objects provided.
    # Note that the 6 coordinates are defined with 3 position in self.r, 3 velocities in self.v.
    
    r = np.array([y[0],y[1],y[2]]) # Create the self position vector
    v = np.array([y[3],y[4],y[5]]) # Create the self velocity vector
    r_dash = v # Trivial ODEs from time derivates
    jup.r = np.array(jup_soln(t)+[0]) # Update Jupiter position in plane z=0 at t
    grav_acc = np.zeros((1,3)) # Make acc vector, note that acceleration will be in (AU s^-2)
    for object in objects:
        dist_vec = object.r-r # Create the vectorial distance
        dist_to_obj = np.linalg.norm(dist_vec) # And the scalar distance
        f_strength = (-4*(np.pi**2)*object.m)/(dist_to_obj**3) # The 'strength' of the force
        grav_acc += f_strength*dist_vec # Acceleration += (-GM/r^3)*_r_, update acceleration
    print(grav_acc) 
    v_dash = grav_acc # Derivative of velocities are the accelerations
    return r_dash.tolist()+v_dash.tolist() # Return the derivatives as concatenated 6-vector

def integrate_gravity(y0, t_points, objects):
    # Function which takes the initial conditions y0 as 6-vector [r1,r2,r3,v1,v2,v3] for 3D space
    # and a series of time points at which to evaluate the solutions it finds (via adaptive
    # stepping) given a landscape of massive bodies in the objects array of Planet class objects.
    # Outputs the y 6-vector at the time points t_points.
    
    y = odeint(grav_derivates,y0,t_points,args=(objects,)) # Calling the integrator function
    return y # Return a list giving the y 6-vector at each time point


