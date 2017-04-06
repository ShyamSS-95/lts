import numpy as np

constants = dict(
                  mass_particle      = 1.0,
                  boltzmann_constant = 1.0,
                )

background = dict(
                  rho      = 1.0,
                  T        = 1.0,
                  vel_bulk = 0
                 )

perturbation = dict(
                    real_part   = 1e-2  ,
                    imag_part   = 0     ,
                    wave_number = 2*np.pi
                   )

size = dict(
            N_vel_x = 1001,
            N_x     = 32,
            vel_max = 10.0
           )

time = dict(
            final_time   = 5.0,
            dt           = 0.01
           )

fields = dict(
             enabled         = 'True',
             charge_particle = 10.0
            )

collisions = dict(
                  enabled = 'True',
                  tau     =  0.005
                 )
