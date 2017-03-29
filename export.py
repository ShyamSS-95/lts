import numpy as np 
import h5py
import lts.calculate as calculate

def density(density_data, time_array):

  h5f = h5py.File('density.h5', 'w')
  h5f.create_dataset('density', data = density_data)
  h5f.create_dataset('time',    data = time_array)
  h5f.close()

def f_distribution(config, final_delta_f, f_background):

  f_dist = np.zeros([N_pos, N_vel])
  for i in range(N_vel):
    f_dist[:, i] = final_delta_f[i] * np.cos(k*x) -\
                   final_delta_f[i + N_vel] * np.sin(k*x)

  f_background = calculate.f
  h5f = h5py.File('f_distribution.h5', 'w')
  h5f.create_dataset('f_dist', data = f_dist)
  h5f.close()