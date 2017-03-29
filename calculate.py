import numpy as np 

class config:
  pass

def set(params):
  config.mass_particle      = params.constants['mass_particle']
  config.boltzmann_constant = params.constants['boltzmann_constant']

  config.rho      = params.background['rho']
  config.T        = params.background['T']
  config.vel_bulk = params.background['vel_bulk']

  config.pert_real_part = params.perturbation['real_part']
  config.pert_imag_part = params.perturbation['imag_part']
  config.wave_number    = params.perturbation['wave_number']

  config.N_vel_x = params.size['N_vel_x']
  config.N_x     = params.size['N_x']
  config.vel_max = params.size['vel_max']

  config.final_time = params.time['final_time']
  config.dt         = params.time['dt']
    
  config.fields_enabled  = params.fields['enabled']
  config.charge_particle = params.fields['charge_particle']

  config.collisions_enabled = params.collisions['enabled']
  config.tau                = params.collisions['tau']

  return config

def f_background(config):
  
  mass_particle      = config.mass_particle
  boltzmann_constant = config.boltzmann_constant

  rho = config.rho
  T   = config.T
  
  vel_max = config.vel_max
  N_vel_x = config.N_vel_x

  vel_x = np.linspace(-vel_max, vel_max, N_vel_x)

  f_background = rho * np.sqrt(mass_particle/(2*np.pi*boltzmann_constant*T)) * \
                 np.exp(-mass_particle*vel_x**2/(2*boltzmann_constant*T))

  return f_background

def diff_f_background(config):

  mass_particle      = config.mass_particle
  boltzmann_constant = config.boltzmann_constant


  T = config.T

  vel_max = config.vel_max
  N_vel_x = config.N_vel_x

  vel_x = np.linspace(-vel_max, vel_max, N_vel_x)

  diff_f_background = f_background(config) *\
                      (-mass_particle*vel_x)/(boltzmann_constant*T)

  return diff_f_background

def time_array(config):

  final_time = config.final_time
  dt         = config.dt

  time_array = np.arange(0, final_time + dt, dt)

  return time_array

def init_delta_fHat(config):

  pert_real_part = config.pert_real_part 
  pert_imag_part = config.pert_imag_part 

  N_vel_x = config.N_vel_x
  vel_max = config.vel_max
  
  vel_x = np.linspace(-vel_max, vel_max, N_vel_x)

  delta_fHat_initial = np.zeros(2 * N_vel_x)

  delta_fHat_initial[:N_vel_x] = pert_real_part * f_background(config) # Real Part
  delta_fHat_initial[N_vel_x:] = pert_imag_part * f_background(config) # Imaginary Part

  return delta_fHat_initial