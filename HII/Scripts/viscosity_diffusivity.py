#!/usr/bin/env python

"""
Calculate the diffusivity and viscosity of solvent in a simulation box using best practices
See : https://github.com/ejmaginn/TransportCheckList/
"""

import argparse
import os
import subprocess
import numpy as np
from genmdp import SimulationMdp
from pymbar.timeseries import detectEquilibration


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-g', '--gro', default='initial.gro', type=str, help='Initial configuration coordinate file')
    parser.add_argument('-lnpt', '--length_npt', default=1, type=float, help='Length of NPT simulation (ns)')
    parser.add_argument('-lnvt', '--length_nvt', default=100, type=float, help='Length of NVT simulation (ns)')
    parser.add_argument('-lnve', '--length_nve', default=100, type=float, help='Length of NVE simulation (ns)')
    parser.add_argument('-T', '--title', default='Generic Molecular Dynamics Simulation', type=str, help='Simulation Title')
    parser.add_argument('-s', '--em_steps', default=5000, type=int, help='Steps to take during energy minimization')
    parser.add_argument('-d', '--dt', default=0.002, type=float, help='time step (ps)')
    parser.add_argument('-p', '--pcoupltype', default='isotropic', type=str, help='Pressure Couple Type')
    parser.add_argument('--temp', default=300, help='Specify temperature at which to run simulation')
    parser.add_argument('--barostat', default='berendsen', type=str, help='pressure coupling scheme to use')
    parser.add_argument('--genvel', default='yes', type=str, help='generate velocities according to a maxwell'
                                                                     'distribution')
    parser.add_argument('--solvent', default='water', help='Name of solvent')
    parser.add_argument('--tau_t', default=1, type=float, help='Temperature coupling time constant')
    parser.add_argument('--tau_p', default=20, type=float, help='Pressure coupling time constant')
    parser.add_argument('-nr', '--nreplicates', default=1, type=int, help='Number of replicate simulations to run')
    parser.add_argument('-f', '--frames', default=50, type=int, help='Number of frames to record. If nstxout, nstvout, '
                        'nstfout, nstenergy are not specified those parameters will be output every frame')
    parser.add_argument('-nx', '--nstxout', type=int, help='Frequency to output coordinates to trajectory file')
    parser.add_argument('-nv', '--nstvout', type=int, help='Frequency to output velocities to trajectory file')
    parser.add_argument('-nf', '--nstfout', type=int, help='Frequency to output forces to trajectory file')
    parser.add_argument('-ne', '--nstenergy', type=int, help='Frequency to output energy to energy file')
    parser.add_argument('-mpi', action="store_true", help='Specify this flag if the system will be run in parallel')
    parser.add_argument('-np', default=4, type=int, help='Number of processes to run in parallel')
    parser.add_argument('-analyze', action="store_true", help='Calculate viscosity and diffusivity only. Do not run'
                                                              'simulations')
    parser.add_argument('-simulate', action="store_true", help='Run simulations only. Do not analyze. This may be good'
                        'to use since the MSD and viscosity calculations will require some manual intervention.')

    args = parser.parse_args()

    return args


class System(object):
    """
    A class to keep track of the system and run various simulations
    """
    def __init__(self, mdp, mpi=False, np=4):
        """
        :param mdp: SimulationMdp object
        :param mpi: Specify if you are running GROMACS in parallel
        :param np: Number of processes. Only has meaning if mpi=True
        """
        self.mdp = mdp
        self.gmx = "gmx"
        if mpi:
            self.gmx = "mpirun -np %d gmx_mpi" % np
        self.top = mdp.top
        self.density = 0
        self.density_equilibration = 0
        self.density_vs_time = []
        self.npt_time = []
        self.nvt_time = []
        self.nve_time = []
        self.replicate_frames = []
        self.mass = self.top.system_mass  # mass of whole system in g/mol

    def energy_minimize(self, configuration, out='em'):
        """
        :param configuration: coordinate file to be minimized
        :param out: name of output file without file extension
        """
        mdp_file = self.mdp.em_mdp_name

        p1 = subprocess.Popen(['gmx', 'grompp', '-f', '%s' % mdp_file, '-p', '%s' % self.mdp.top.name, '-o',
                              '%s' % out, '-c', '%s' % configuration])  # will make ensemble.trr, ensemble.gro etc.
        p1.wait()
        run_string = "%s mdrun -v -deffnm %s" % (self.gmx, out)
        p2 = subprocess.Popen(run_string, shell=True)
        # p2 = subprocess.Popen(['%s' % self.gmx, 'mdrun', '-v', '-deffnm', '%s' % out])
        p2.wait()

    def run_simulation(self, ensemble, configuration):
        """
        :param ensemble: thermodynamic ensemble to run simulation in
        :param configuration: initial configuration coordinate file
        """

        mdp_file = None  # dummy so pycharm stops yelling at me
        if ensemble == 'nvt':
            mdp_file = self.mdp.nvt_mdp_name
        elif ensemble == 'npt':
            mdp_file = self.mdp.npt_mdp_name
        elif ensemble == 'nve':
            mdp_file = self.mdp.nve_mdp_name
        else:
            print('Please select a valid ensemble to run your simulation')
            exit()

        p1 = subprocess.Popen(['gmx', 'grompp', '-f', '%s' % mdp_file, '-p', '%s' % self.mdp.top.name, '-o',
                              '%s' % ensemble, '-c', '%s' % configuration])  # will make ensemble.trr, ensemble.gro etc.
        p1.wait()

        run_string = "%s mdrun -v -deffnm %s" % (self.gmx, ensemble)
        p2 = subprocess.Popen(run_string, shell=True)
        # p2 = subprocess.Popen(['%s' % self.gmx, 'mdrun', '-v', '-deffnm', '%s' % ensemble])
        p2.wait()

    def average_density(self, ensemble):
        """
        :param ensemble: Thermodynamic ensemble that simulation was run in. Will dictate names of files
        :return: average density after equilibration
        """

        # get density vs. time using gmx energy
        ps = subprocess.Popen(('echo', 'Density'), stdout=subprocess.PIPE)
        ps.wait()
        p = subprocess.Popen(('gmx', 'energy', '-f', '%s.edr' % ensemble), stdin=ps.stdout,
                             stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        p.wait()

        # read output file (default is energy.xvg. This will always be the output file so it's hardcoded for now)
        with open('energy.xvg', 'r') as f:
            a = []
            for line in f:
                a.append(line)

        data_start = 0
        while a[data_start].count('Density') == 0:
            data_start += 1

        data_start += 1
        d = np.zeros([len(a) - data_start])  # density
        t = np.zeros_like(d)  # time

        for i in range(data_start, len(a)):
            line_data = a[i].split()
            t[i - data_start] = line_data[0]
            d[i - data_start] = line_data[1]

        # detect when data decorrelates from itself
        self.density_equilibration = detectEquilibration(d[10:])[0]  # discard first few data points to avoid pymbar bug.
        print('Density equilibrated after %d frames' % self.density_equilibration)
        self.density = np.mean(d[self.density_equilibration:])  # GROMACS outputs in kg / m^3
        self.density_vs_time = d  # save density vs. time data

        if ensemble == 'nvt':
            self.nvt_time = t
        elif ensemble == 'npt':
            self.npt_time = t
        elif ensemble == 'nve':
            self.nve_time = t

        return self.density

    def npt_replicates(self, nreplicates):
        """
        :param nreplicates: Number of simulation replicates that will be run
        """

        choices = [i for i in range(len(self.npt_time))]
        choices = choices[self.density_equilibration:]  # choose configuration from any time point after equilibration
        self.replicate_frames = np.random.choice(choices, size=nreplicates, replace=False)

    def generate_replicate_gro(self, n, ensemble='npt', out='replicate.gro'):
        """
        :param n: replicate number you want (refers to index in self.replicate_frames) (int)
        :param ensemble: thermodynamic ensemble from which to pull replicate (str)
        :param out: name of output file (str)
        """

        t = self.npt_time[self.replicate_frames[n]]

        ps = subprocess.Popen(('echo', 'System'), stdout=subprocess.PIPE)  # want coordinates of full system
        ps.wait()
        p1 = subprocess.Popen(('gmx', 'trjconv', '-f', '%s.trr' % ensemble, '-s', '%s.tpr' % ensemble,
                              '-o', '%s' % out, '-b', '%s' % t, '-dump', '%s' % t), stdin=ps.stdout,
                             stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # suppress output
        p1.wait()

        NA = 6.022 * 10 ** 23  # avogradros number
        m3_2_nm3 = 1 * 10 ** 27  # convert meters cubed to nanometers cubed
        box_volume = m3_2_nm3 * self.mass / (1000 * NA * self.density)
        cubic_box_length = box_volume ** (1/3)

        # make new box with same density as average of NPT simulation
        p2 = subprocess.Popen(('gmx', 'editconf', '-f', '%s' % out, '-bt', 'cubic', '-box', '%s' % cubic_box_length,
                               '%s' % cubic_box_length,  '%s' % cubic_box_length, '-o', '%s' % out))
        p2.wait()


if __name__ == "__main__":

    args = initialize()

    if not args.analyze:

        # get output frequencies (important for controlling size of trajectory)
        if args.nstxout:
            nstxout = args.nstxout
        else:
            nstxout = int(args.length_npt*1000 / (args.dt * args.frames))

        if args.nstvout:
            nstvout = args.nstvout
        else:
            nstvout = int(args.length_npt*1000 / (args.dt * args.frames))

        if args.nstfout:
            nstfout = args.nstfout
        else:
            nstfout = int(args.length_npt*1000 / (args.dt * args.frames))

        if args.nstenergy:
            nstenergy = args.nstenergy
        else:
            nstenergy = int(args.length_npt*1000 / (args.dt * args.frames))  # output frequency

        mdp = SimulationMdp(args.gro, title=args.title, T=args.temp, em_steps=args.em_steps,
                            time_step=args.dt, length=args.length_npt*1000, p_coupling=args.pcoupltype,
                            barostat=args.barostat, genvel=args.genvel, tau_p=args.tau_p,
                            tau_t=args.tau_t, nstxout=nstxout, nstenergy=nstenergy, nstvout=nstvout, nstfout=nstfout)

        top = mdp.top  # SimulationMdp calls to SystemTopology. Redefined simply as top for readability
        top.write_top()  # write system topology

        mdp.write_em_mdp()  # write energy minimization .mdp

        mdp.nstenergy = int(args.length_npt * 1000 / (args.dt * 1000))  # let's get 1000 energy outputs for density calculation
        mdp.write_npt_mdp()  # write .mdp for npt simulation  -- PARRINELLO-RAHMAN?? TEST STABILITY

        mdp.length = args.length_nvt * 1000  # convert to picoseconds
        mdp.nstenergy = int(np.ceil(5 / (args.dt*1000)))  # save every 6 femptoseconds (3 timesteps!)
        mdp.write_nvt_mdp()  # write .mdp for nvt simulation

        mdp.length = args.length_nve * 1000  # convert to picoseconds
        mdp.nstenergy = int(np.ceil(5 / (args.dt*1000)))  # save every 6 femptoseconds (3 timesteps!)
        mdp.write_nve_mdp()  # write .mdp for nve simulation

        if args.simulate:

            sys = System(mdp, mpi=args.mpi, np=args.np)  # system object to keep track of all relevant data

            # Run NPT equilibraton simulation
            sys.run_simulation('npt', args.gro)
            average_density = sys.average_density('npt')
            sys.npt_replicates(args.nreplicates)  # decide which frames to use as starting points for replicate simulations
            sys.generate_replicate_gro(0)  # create .gro with average density from npt equilibration

            #sys.energy_minimize('replicate.gro', out='replicate_em')
            sys.run_simulation('nvt', 'replicate.gro')
            #sys.run_simulation('nve', 'nvt.gro')  # run nve simulation using last frame of nvt simulation

    # Now calculate diffusivity and viscosity
    # if not args.simulate:
    #
    #     V = viscosity.Viscosity(args.edr, args.gro)
    #     V.plot_all()
