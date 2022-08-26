from Integrators import euler_osc_solver, verlet_osc_solver, verlet_ballistic_solver
import matplotlib.pyplot as plt
from numpy import pi
import numpy as np

if __name__ == "__main__":

    # Running the harmonic oscillator - commented
    '''
    # Parameters for the integrator
    time_params = [ 0.0 , 10*pi , 1000 ]
    pos_i = 1.0
    vel_i = 0.0
    spring_params = [ 1.0 , 1.0 , 0.0 ]

    time , pos , vel = euler_osc_solver( time_params , pos_i , vel_i , spring_params )
    time , pos_v , vel_v = verlet_osc_solver( time_params , pos_i , vel_i , spring_params )

    # Real solution to compare with
    pos_r = np.cos( time )

    # 2D Plot example of integrated vs. real solution
    fig, ax = plt.subplots()
    ax.plot( time , pos , color = "green" , linestyle = "solid" , label = r"result Euler" )
    ax.plot( time , pos_r , color = "red" , linestyle = "solid" , label = r"real position" )
    ax.plot( time , pos_v , color = "blue" , linestyle = "solid" , label = r"result Verlet" )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"Position in [m]" )

    ax.legend( loc = "upper right" )

    #ax.set_ylim( 0.0 , pi/2.0 )
    plt.grid( True )

    #fig.savefig( "fig_name.pdf" , format = "pdf" )
    
    plt.show()
    '''
    # Running the ballistic propagator
    # Parameters for the integrator
    time_params = [ 0.0 , 1.8 , 1000 ]
    vel_tot = 10.0 # [m/s]
    alpha = 60.0 # [deg]
    ball_params = 0.0

    time , pos , vel = verlet_ballistic_solver( time_params , vel_tot , alpha , ball_params )

    # Transpose positions to get time arrays separately
    x_plot = np.transpose( pos )[ 0 ]
    y_plot = np.transpose( pos )[ 1 ]

    # 2D Plot example of integrated vs. real solution
    fig, ax = plt.subplots()
    ax.plot( x_plot , y_plot , color = "green" , linestyle = "solid" , label = r"Trajectory" )

    ax.set_xlabel( r"X position [m]" )
    ax.set_ylabel( r"Y position [m]" )

    ax.legend( loc = "upper right" )
    plt.grid( True )

    #fig.savefig( "fig_name.pdf" , format = "pdf" )
    
    plt.show()