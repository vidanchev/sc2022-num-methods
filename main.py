from Integrators import euler_osc_solver, verlet_osc_solver
import matplotlib.pyplot as plt
from numpy import pi
import numpy as np

if __name__ == "__main__":

    # Parameters for the integrator
    time_params = [ 0.0 , 10*pi , 1000 ]
    pos_i = 1.0
    vel_i = 0.0
    spring_params = [ 1.0 , 1.0 ]

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