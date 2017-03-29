import lts.calculate as calculate
import lts.evolve as evolve
import lts.plot as plot
import params

config = calculate.set(params)

delta_f_hat_initial = calculate.init_delta_fHat(config)
time_array          = calculate.time_array(config)

rho, f_final = evolve.solve(config, delta_f_hat_initial, time_array)

plot.density(rho, time_array)