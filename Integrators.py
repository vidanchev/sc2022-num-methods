import numpy as np

g_const = 9.8 # Gravitational constant in [m/s^2]
beta = 1.5 # Ballistic coefficient [1/s]

# Right Hand Side function for Oscillator
# Inputs:
# - spring mass k_spring [N/m]
# - mass [kg]
# - beta [N*s/m] friction coefficient
# - pos [m]: position (1D)
# - vel [m/s]: velocity (1D)
# Outputs:
# - force [N]
def rhs_oscillator( k_spring , mass , beta , pos , vel ):

    # Oscillator with friction -> F = - k*x - beta*v
    force =  - k_spring*pos/mass - beta*vel

    return force

# Euler solver
# Inputs:
# - time_params[ 3 ]: [ t_i , t_f , npoints ]
# -- t_i [sec]: initial time
# -- t_f [sec]: final time
# -- npoints [-]: number of points
# - pos_i [m]: initial position (1D)
# - vel_i [m/s]: initial velocity (1D)
# - spring_params[ 2 ]: [ k_spring , mass , beta ]
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
    beta = spring_params[ 2 ]
    
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
        vel[ i + 1 ] = vel[ i ] + rhs_oscillator( k_spring , mass , beta , pos[ i ] , vel[ i ] )*dt
    
    return time , pos , vel

# Verlet solver
# Inputs:
# - time_params[ 3 ]: [ t_i , t_f , npoints ]
# -- t_i [sec]: initial time
# -- t_f [sec]: final time
# -- npoints [-]: number of points
# - pos_i [m]: initial position (1D)
# - vel_i [m/s]: initial velocity (1D)
# - spring_params[ 2 ]: [ k_spring , mass , beta ]
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
    beta = spring_params[ 2 ]
    
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
        v_half = vel[ i ] + rhs_oscillator( k_spring , mass , beta , pos[ i ] , vel[ i ] )*dt/2.0
        pos[ i + 1 ] = pos[ i ] + v_half*dt
        vel[ i + 1 ] = v_half + rhs_oscillator( k_spring , mass , beta , pos[ i + 1 ] , v_half )*dt/2.0
    
    return time , pos , vel

# Inputs:
# - pos[ 2 ] [m]: position array [ x , y ]
# - vel[ 2 ] [m/s]: velocity array [ x , y ]
# Outputs:
# - rhs_pos[ 2 ]
# - rhs_vel[ 2 ]
def rhs_ballistic( pos , vel ):
    
    # right-hand side of position and velocity
    rhs_pos = [ vel[ 0 ] , vel[ 1 ] ] # this is the right-hand-side for position [ 0 ] and [ 1 ] are "x" and "y"
    rhs_vel = [ - beta*vel[ 0 ] , - g_const - beta*vel[ 1 ] ] # this is the right-hand-side for velocity [ 0 ] and [ 1 ] are "x" and "y"

    return rhs_pos , rhs_vel

# Verlet solver for ballistic trajectory
# Inputs:
# - time_params[ 3 ]: [ t_i , t_f , npoints ]
# -- t_i [sec]: initial time
# -- t_f [sec]: final time
# -- npoints [-]: number of points
# - vel_tot [m/s]: total velocity
# - alpha [deg]: angle of "hit"
# Outputs:
# - time[ npoints ] [sec]: time array
# - pos[ npoints ][ 2 ] [m]: position array [ x , y ]
# - vel[ npoints ][ 2 ] [m/s]: velocity array [ x , y ]
def verlet_ballistic_solver( time_params , vel_tot , alpha , ball_params ):

    # Unpackage the parameters
    t_i = time_params[ 0 ]
    t_f = time_params[ 1 ]
    npoints = int( time_params[ 2 ] )

    # Find my timestep by dividing by number of intervals
    dt = ( t_f - t_i )/( npoints - 1 )

    # Make a time array starting from t_i to t_f with npoints (number of points)
    time = np.linspace( t_i , t_f , npoints )
    # make position and velocity arrays which are empty -> we have to populate them
    # arrays are [ npoints ][ 2 ] -> [ i ][ 0 ] is x[ i ] and [ i ][ 1 ] is y[ i ]
    pos = np.zeros( ( npoints , 2 ) )
    vel = np.zeros( ( npoints , 2 ) )

    pos[ 0 ] = [ 0.0 , 0.0 ]
    vel[ 0 ] = [ vel_tot*np.cos( alpha*np.pi/180.0 ) , vel_tot*np.sin( alpha*np.pi/180.0 ) ]

    # Main integration loop
    # Find pos[ i + 1 ][ j ] and v[ i + 1 ][ j ] knowing the current ones ( pos[ i ] and vel[ i ] )
    for i in range( 0 , npoints - 1 ):

        # This is the better way to do the integrator -> looping over indices instead of writing everything for each dimension
        # Find the half-velocity step
        rhs_pos, rhs_vel = rhs_ballistic( pos[ i ] , vel[ i ] )
        v_half = np.zeros( 2 )
        for j in range( 0 , 2 ):
            v_half[ j ] = vel[ i ][ j ] + rhs_vel[ j ]*dt/2.0
        
        # Find the next position step based on half-velocity
        rhs_pos, rhs_vel = rhs_ballistic( pos[ i ] , v_half )
        for j in range( 0 , 2 ):
            pos[ i + 1 ][ j ] = pos[ i ][ j ] + rhs_pos[ j ]*dt
        
        # Find the next velocity step
        rhs_pos, rhs_vel = rhs_ballistic( pos[ i + 1 ] , v_half )
        for j in range( 0 , 2 ):
            vel[ i + 1 ][ j ] = v_half[ j ] + rhs_vel[ j ]*dt/2.0

        # This is the "amateur" way to do the integrator -> adding manually all steps
        '''
        # Find the half-velocity step
        v_half_x = vel[ i ][ 0 ] + ( 0.0 )*dt/2.0
        v_half_y = vel[ i ][ 1 ] + ( - g_const )*dt/2.0

        # Find the next position step based on half-velocity
        pos[ i + 1 ][ 0 ] = pos[ i ][ 0 ] + v_half_x*dt
        pos[ i + 1 ][ 1 ] = pos[ i ][ 1 ] + v_half_y*dt

        # Find the next velocity step
        vel[ i + 1 ][ 0 ] = v_half_x + ( 0.0 )*dt/2.0
        vel[ i + 1 ][ 1 ] = v_half_y + ( - g_const )*dt/2.0
        '''

        # simple way of stopping the integration -> if we "hit the ground" ( y = 0 ) -> we keep y = 0 (we roll on the surface)
        if pos[ i + 1 ][ 1 ] < 0.0:
            pos[ i + 1 ][ 1 ] = 0.0

    return time , pos , vel