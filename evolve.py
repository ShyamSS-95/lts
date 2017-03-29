import numpy as np
from scipy.integrate import odeint
import lts.calculate as calculate

def linearized_collision_operator(config, f):

  mass_particle      = config.mass_particle
  boltzmann_constant = config.boltzmann_constant

  T   = config.T
  rho = config.rho

  vel_max = config.vel_max
  N_vel_x = config.N_vel_x

  tau   = config.tau

  vel_x = np.linspace(-vel_max, vel_max, N_vel_x)
  dv    = vel_x[1] - vel_x[0]

  delta_T   = np.sum(f * (vel_x**2 - T)) * dv/rho
  delta_rho = np.sum(f) * dv
  delta_v   = np.sum(f * vel_x) * dv/rho
  
  expr_term_1 = np.sqrt(2 * mass_particle**3) * delta_T * rho * vel_x**2
  expr_term_2 = 2 * np.sqrt(2 * mass_particle) * boltzmann_constant * delta_rho * T**2
  expr_term_3 = 2 * np.sqrt(2 * mass_particle**3) * rho * delta_v * vel_x * T
  expr_term_4 = - np.sqrt(2 * mass_particle) * boltzmann_constant * delta_T * rho * T
  
  C_f = (((expr_term_1 + expr_term_2 + expr_term_3 + expr_term_4)*\
         np.exp(-mass_particle * vel_x**2/(2 * boltzmann_constant * T))/\
         (4 * np.sqrt(np.pi * T**5 * boltzmann_constant**3)) - f
         )/tau
        )
  
  return C_f

def diff_delta_f(Y, t, config):

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

  f_r = Y[:vel_x.size]
  f_i = Y[vel_x.size:]

  int_Df_i = np.sum(f_i) * dv
  int_Df_r = np.sum(f_r) * dv

  if(fields_enabled!="True"):

    fields_term_r = 0
    fields_term_i = 0

  else:

    int_Df_i = np.sum(f_i) * dv
    int_Df_r = np.sum(f_r) * dv

    diff_f_background = calculate.diff_f_background(config)

    fields_term_r = - charge_particle**2 *(int_Df_i * diff_f_background/wave_number)
    fields_term_i = + charge_particle**2 *(int_Df_r * diff_f_background/wave_number)

  if(collisions_enabled!="True"):

    C_f_r = 0
    C_f_i = 0

  else:

    f     = f_r + 1j * f_i
    C_f   = linearized_collision_operator(config, f)
    C_f_r = C_f.real
    C_f_i = C_f.imag

  dYdt =np.concatenate([(wave_number * vel_x * f_i)   - fields_term_r + C_f_r,\
                         -(wave_number * vel_x * f_r) + fields_term_i + C_f_i\
                       ], axis = 0)
  
  return dYdt

def solve(config, delta_f_initial, time_array):
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
    
    rho_r = np.sum(new_delta_f[:N_vel_x])*dv
    rho_i = np.sum(new_delta_f[N_vel_x:])*dv

    density_data[time_index] = np.max(rho_r*np.cos(k*x) - rho_i*np.sin(k*x))
    old_delta_f              = new_delta_f.copy()

  return(density_data, new_delta_f)