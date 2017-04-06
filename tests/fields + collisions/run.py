import lts.calculate as calculate
import lts.evolve as evolve
# import lts.plot as plot
import pylab as pl

import tau_1e2_
import tau_1e3_
import tau_5e2_
import tau_5e3_
import tau_inf_

config = calculate.set(tau_1e2_)

delta_f_hat_initial = calculate.init_delta_fHat(config)
time_array          = calculate.time_array(config)

rho_1, f_final = evolve.solve(config, delta_f_hat_initial, time_array)

config = calculate.set(tau_5e2_)

delta_f_hat_initial = calculate.init_delta_fHat(config)
time_array          = calculate.time_array(config)

rho_2, f_final = evolve.solve(config, delta_f_hat_initial, time_array)


config = calculate.set(tau_5e3_)

delta_f_hat_initial = calculate.init_delta_fHat(config)
time_array          = calculate.time_array(config)

rho_3, f_final = evolve.solve(config, delta_f_hat_initial, time_array)


config = calculate.set(tau_1e3_)

delta_f_hat_initial = calculate.init_delta_fHat(config)
time_array          = calculate.time_array(config)

rho_4, f_final = evolve.solve(config, delta_f_hat_initial, time_array)


config = calculate.set(tau_inf_)

delta_f_hat_initial = calculate.init_delta_fHat(config)
time_array          = calculate.time_array(config)

rho_5, f_final = evolve.solve(config, delta_f_hat_initial, time_array)

pl.plot(time_array, rho_1, label = r'$\tau=0.01$')
pl.plot(time_array, rho_2, label = r'$\tau=0.05$')
pl.plot(time_array, rho_3, label = r'$\tau=0.005$')
pl.plot(time_array, rho_4, label = r'$\tau=0.001$')
pl.plot(time_array, rho_5, label = r'$\tau=\infty$')
pl.legend()
pl.show()
