import numpy as np

# Right Hand Side function for Oscillator
# Inputs:
# - spring mass k_spring [N/m]
# - mass [kg]
# - pos [m]: position (1D)
# - vel [m/s]: velocity (1D)
# Outputs:
# - force [N]
def rhs_oscillator( k_spring , mass , pos , vel ):

    force =  - k_spring*pos/mass

    return force

# Euler solver
# Inputs:
# - time_params[ 3 ]: [ t_i , t_f , npoints ]
# -- t_i [sec]: initial time
# -- t_f [sec]: final time
# -- npoints [-]: number of points
# - pos_i [m]: initial position (1D)
# - vel_i [m/s]: initial velocity (1D)
# - spring_params[ 2 ]: [ k_spring , mass ]
# Outputs:
# - time[ npoints ] [sec]: time array
# - pos[ npoints ] [m]: position array
# - vel[ npoints ] [m/s]: velocity array
def euler_osc_solver( time_params , pos_i , vel_i , spring_params ):

    # Unpackage the parameters
    t_i = time_params[ 0 ]
    t_f = time_params[ 1 ]
    npoints = int( time_params[ 2 ] )
    k_spring = spring_params[ 0 ]
    mass = spring_params[ 1 ]
    
    # Find my timestep by dividing by number of intervals
    dt = ( t_f - t_i )/( npoints - 1 )

    # Make a time array starting from t_i to t_f with npoints (number of points)
    time = np.linspace( t_i , t_f , npoints )
    # make position and velocity arrays which are empty -> we have to populate them
    pos = np.zeros( npoints )
    vel = np.zeros( npoints )

    # Populate with initial position and velocity
    pos[ 0 ] = pos_i
    vel[ 0 ] = vel_i

    # Main integration loop
    # Find x[ i + 1 ] and v[ i + 1 ] knowing the current ones ( x[ i ] and v[ i ] )
    for i in range( 0 , npoints - 1 ):
        pos[ i + 1 ] = pos[ i ] + vel[ i ]*dt
        vel[ i + 1 ] = vel[ i ] + rhs_oscillator( k_spring , mass , pos[ i ] , vel[ i ] )*dt
    
    return time , pos , vel

# Verlet solver
# Inputs:
# - time_params[ 3 ]: [ t_i , t_f , npoints ]
# -- t_i [sec]: initial time
# -- t_f [sec]: final time
# -- npoints [-]: number of points
# - pos_i [m]: initial position (1D)
# - vel_i [m/s]: initial velocity (1D)
# - spring_params[ 2 ]: [ k_spring , mass ]
# Outputs:
# - time[ npoints ] [sec]: time array
# - pos[ npoints ] [m]: position array
# - vel[ npoints ] [m/s]: velocity array
def verlet_osc_solver( time_params , pos_i , vel_i , spring_params ):

    # Unpackage the parameters
    t_i = time_params[ 0 ]
    t_f = time_params[ 1 ]
    npoints = int( time_params[ 2 ] )
    k_spring = spring_params[ 0 ]
    mass = spring_params[ 1 ]
    
    # Find my timestep by dividing by number of intervals
    dt = ( t_f - t_i )/( npoints - 1 )

    # Make a time array starting from t_i to t_f with npoints (number of points)
    time = np.linspace( t_i , t_f , npoints )
    # make position and velocity arrays which are empty -> we have to populate them
    pos = np.zeros( npoints )
    vel = np.zeros( npoints )

    # Populate with initial position and velocity
    pos[ 0 ] = pos_i
    vel[ 0 ] = vel_i

    # Main integration loop
    # Find x[ i + 1 ] and v[ i + 1 ] knowing the current ones ( x[ i ] and v[ i ] )
    for i in range( 0 , npoints - 1 ):
        v_half = vel[ i ] + rhs_oscillator( k_spring , mass , pos[ i ] , vel[ i ] )*dt/2.0
        pos[ i + 1 ] = pos[ i ] + v_half*dt
        vel[ i + 1 ] = v_half + rhs_oscillator( k_spring , mass , pos[ i + 1 ] , vel[ i + 1 ] )*dt/2.0
    
    return time , pos , vel


def prop_func( ):

    print( "This will be a propagator function" )