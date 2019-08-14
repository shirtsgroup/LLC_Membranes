#!/usr/bin/env python

from LLC_Membranes.analysis.identify_states import States, Chain
from LLC_Membranes.llclib import file_rw
from scipy.stats import levy_stable

solute = 'MET'
directory = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'

states = file_rw.load_object('%s/%s/10wt/states.pl' % (directory, solute))

# save emissions plot
states.plot_emissions(save=True, savename='emissions_%s' % solute)

chains = Chain(states.transition_matrix, states.fit_params, emission_function=levy_stable)
chains.generate_realizations(24, 2000, bound=states.box_length)
chains.calculate_msd()
chains.plot_msd(save=True, savename='%s_msd.pdf' % solute)

