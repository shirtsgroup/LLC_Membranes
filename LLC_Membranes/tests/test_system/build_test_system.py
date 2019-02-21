#!/usr/bin/env python

from LLC_Membranes.setup import solvation_equilibration as equil
import subprocess

args = equil.initialize().parse_args()  # use default parameters. Change a few below
args.monomers_per_column = 6

sys = equil.System(args)
sys.r = 0.6 # database entry which give ~ 10 wt % water
sys.equilibrate()
sys.calculate_pore_water()
sys.write_final_pore_configuration()
sys.place_water_tails()
