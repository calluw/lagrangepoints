# File which uses an adaptive LSODA method in order to investigate the orbits of
# negligible masses around the Lagrange points of large mass planet-star systems, to 
# work out the scaling relationship of the stable orbit region with the planetary mass.
# Units used in the problem are length unit 1AU, time unit 1yr, mass unit 1solarmass
# In this file, the frame of reference is the frame rotating about the barycentre of
# the Sun-Jupiter system at the rotational speed of Jupiter, hence includes inertial
# forces

## Imports
from main import * # Take existing objects and barycentre calculation.

def grav_derivatives_rot(y, t, objects):
    # Function which calculates the derivatives of the planet's position and velocity
    # at a given time t, using the gravitational forces from the objects provided.
    # Note that the 6 coordinates are defined with 3 position in self.r, 3 velocities in self.v.
    # Includes the centrifugal force from the rotating frame
    
    global jup, dist # Import Jupiter object to work out speed of rotation of frame
    jup_r = np.linalg.norm(jup.r0)
    r = np.array([y[0],y[1],y[2]]) # Create the self position vector
    v = np.array([y[3],y[4],y[5]]) # Create the self velocity vector
    r_dash = v # Trivial ODEs from time derivates
    grav_acc = np.zeros(3) # Make acc vector, note that acceleration will be in (AU s^-2)
    for object in objects:
        dist_vec = r-object.r0 # Create the vectorial distance, using the fixed Sun/Jupiter
        dist_to_obj = np.linalg.norm(dist_vec) # And the scalar distance
        f_strength = (-4*(np.pi**2)*object.m)/(dist_to_obj**3) # The 'strength' of the force
        grav_acc += f_strength*dist_vec # Acceleration += (-GM/r^3)*_r_, update acceleration
    trojan_r = np.linalg.norm(r)
    trojan_angle = np.arctan2(r[1],r[0]) # Get the angle of the Trojan for centrifugal force
    frame_angular_speed = 2*np.pi*np.sqrt(1/(jup_r*(dist**2)))
    centrifugal_acc_strength =  trojan_r*(frame_angular_speed**2) # Centrifugal is r*(omega)^2
    centrifugal_acc = centrifugal_acc_strength * np.array([np.cos(trojan_angle), np.sin(
        trojan_angle), 0]) # And it always points radially out from the rotation point
    coriolis_acc = -2*np.cross([0,0,frame_angular_speed],v) # Coriolis is -2*omega^v
    v_dash = grav_acc + centrifugal_acc + coriolis_acc # Derivative of velocities are accelerations
    return r_dash.tolist()+v_dash.tolist() # Return the derivatives as concatenated 6-vector

def integrate_gravity_rot(y0, t_points, objects):
    # Function which takes the initial conditions y0 as 6-vector [r1,r2,r3,v1,v2,v3] for 3D space
    # and a series of time points at which to evaluate the solutions it finds (via adaptive
    # stepping) given a landscape of massive bodies in the objects array of Planet class objects.
    # Outputs the y 6-vector at the time points t_points.
    
    y = odeint(grav_derivatives_rot,y0,t_points,args=(objects,)) # Calling the integrator function
    return y # Return a list giving the y 6-vector at each time point

def trojan_stability_measure_rot(t_points, soln, sign):
    # Function which takes the end point of the asteroid at the final time point
    # provided, and then calculates the exact Jupiter point at that time, and then
    # the theoretical exact position of the Lagrance point being used (given that
    # it was at sgn(sign)*pi/3 from Jupiter-like), then returns the difference as
    # a fraction of the Jupiter radius. 
    global jup, sun, dist # Read in Jupiter, Sun objects and Sun-Jupiter distance
    trojan_soln = np.array(soln[-1][0:3]) # Gives the calculated Trojan position at final t_points
    jup_r = np.linalg.norm(jup.r0)
    trojan_angle = sign*(np.pi/3) # Generate angle planet would be at
    theor_trojan_soln = sun.r0 + np.array([dist*np.cos(trojan_angle),dist*np.sin(trojan_angle),0])
    trojan_diff_dist = np.linalg.norm(trojan_soln-theor_trojan_soln)
    trojan_frac_diff = trojan_diff_dist/jup_r
    return trojan_frac_diff

def perturbed_initial_condition_rot(dr,dv,angle):
    # Function which produces the initial condition 6-vector for a planet
    # y0 = [x,y,z,vx,vy,vz] based on its radian  angle from Jupiter (negative is counted
    # as clockwise from Jupiter) and a perturbation (dr, dv) from that point.
    global jup, sun # Read in the celestial objects
    global dist # Read in the solar-Jupiter distance
    jup_r = np.linalg.norm(jup.r) # Get the radius of Jupiter-like
    jup_angle = np.arcsin(jup.r0[1]/jup_r) # Initial angle of Jupiter-like
    new_angle = jup_angle + angle # Angle wrt the Sun
    planet_x = dist*np.cos(new_angle)+sun.r0[0] # Calculate the angle from Sun position!
    planet_y = dist*np.sin(new_angle)+sun.r0[1]
    planet_pos = [planet_x+dr[0], planet_y+dr[1], dr[2]] # Unperturbed pos is in plane z=0
    # Now return 6-vector:
    return planet_pos+dv

