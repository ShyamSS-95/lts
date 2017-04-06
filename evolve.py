import numpy as np
from scipy.integrate import odeint
import lts.initialize as initialize
from lts.collision_operators import BGK_collision_operator

def ddelta_fHat_dt(Y, t, config):

  mass_particle      = config.mass_particle
  boltzmann_constant = config.boltzmann_constant

  rho = config.rho
  T   = config.T
  
  vel_max = config.vel_max
  N_vel_x = config.N_vel_x
  
  vel_x = np.linspace(-vel_max, vel_max, N_vel_x)
  dv    = vel_x[1] - vel_x[0]

  wave_number    = config.wave_number   
  
  fields_enabled  = config.fields_enabled
  charge_particle = config.charge_particle

  collisions_enabled = config.collisions_enabled 
  tau                = config.tau

  delta_f_hat_r = Y[:vel_x.size]
  delta_f_hat_i = Y[vel_x.size:]

  delta_f_hat   = delta_f_hat_real + 1j*delta_f_hat_imag

  delta_rho_hat = np.sum(delta_f_hat) * dv

  if(fields_enabled!="True"):

    fields_term_r = 0
    fields_term_i = 0

  else:

    dfdv_background = calculate.diff_f_background(config)

    delta_E_hat =   charge_particle # MAKE SURE THIS IS CORRECT 
                  * (delta_rho) / (1j * kx)

    fields_term =  charge_particle / mass_particle 
                  * delta_E_hat * dfdv_background_

  if(collisions_enabled!="True"):

    C_f_r = 0
    C_f_i = 0

  else:

    f     = delta_f_hat_r + 1j * delta_f_hat_i
    C_f   = linearized_collision_operator(config, f)
    C_f_r = C_f.real
    C_f_i = C_f.imag

  dYdt =np.concatenate([(wave_number * vel_x * delta_f_hat_i)   -
                         fields_term.real + C_f_r,\
                         -(wave_number * vel_x * delta_f_hat_r) +
                         fields_term.imag + C_f_i\
                       ], axis = 0)
  
  return dYdt

def time_integration(config, delta_f_initial, time_array):
  
  vel_max = config.vel_max
  N_vel_x = config.N_vel_x
  N_x     = config.N_x

  x       = np.linspace(0, 1, N_x)
  
  k       = config.wave_number 
  vel_x   = np.linspace(-vel_max, vel_max, N_vel_x)
  dv      = vel_x[1] - vel_x[0]  

  density_data = np.zeros(time_array.size)

  for time_index, t0 in enumerate(time_array):
    t0 = time_array[time_index]
    if (time_index == time_array.size - 1):
        break
    t1 = time_array[time_index + 1]
    t = [t0, t1]

    if(time_index == 0):
        initial_conditions_delta_f = delta_f_initial.copy()
        
    else:
        initial_conditions_delta_f = old_delta_f.copy()
        
    # Integrating delta f
    new_delta_f = odeint(diff_delta_f, initial_conditions_delta_f, t, args = (config,),
                         rtol = 1e-20, atol = 1e-15
                        )[1]
    
    delta_rho_hat_r = np.sum(new_delta_f[:N_vel_x])*dv
    delta_rho_hat_i = np.sum(new_delta_f[N_vel_x:])*dv

    density_data[time_index] = np.max(delta_rho_hat_r * np.cos(k*x) - delta_rho_hat_i * np.sin(k*x))
    old_delta_f              = new_delta_f.copy()

  return(density_data, new_delta_f)
