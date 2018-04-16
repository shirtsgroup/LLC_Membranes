#!/usr/bin/env python

import argparse
from gentop import SystemTopology


def initialize():

    parser = argparse.ArgumentParser(description='Generate .mdp file for simulation of choice')

    parser.add_argument('-T', '--title', default='MD Simulation', type=str, help='Simulation Title')
    parser.add_argument('-g', '--gro', default='initial.gro', type=str, help='coordinate file of system to be simulated')
    parser.add_argument('-t', '--itp', default='dipole.itp', type=str, help='Name of .itp describing monomers')
    parser.add_argument('-s', '--em_steps', default=5000, type=int, help='Steps to take during energy minimization')
    parser.add_argument('-e', '--ensemble', default='npt', type=str, help='Thermodynamic ensemble to put system in')
    parser.add_argument('-d', '--dt', default=0.002, type=float, help='time step (ps)')
    parser.add_argument('-l', '--length', default=1000, type=float, help='simulation length (ps)')
    parser.add_argument('-f', '--frames', default=500, type=int, help='number of frames')
    parser.add_argument('-p', '--pcoupltype', default='semiisotropic', type=str, help='Pressure Couple Type')
    parser.add_argument('--restraints', help='If restraints are on, another mdp option needs to be turned on, so specify '
                                             'this flag', action="store_true")
    parser.add_argument('-x', '--xlink', help='Turn this to "on" if the the system is crosslinked', action="store_true")
    parser.add_argument('-S', '--solvate', help='Specify this if the system has water so an extra line can be added to the '
                                                'topology', action="store_true")
    parser.add_argument('--temp', default=300, help='Specify temperature at which to run simulation')
    parser.add_argument('--mdp', action="store_true", help='Only the .mdp will be written if this option is specified')
    parser.add_argument('--barostat', default='berendsen', type=str, help='pressure coupling scheme to use')
    parser.add_argument('--genvel', default='yes', type=str, help='generate velocities according to a maxwell'
                                                                     'distribution')
    parser.add_argument('--bcc', action="store_true", help='Generate input files using bicontinuous cubic files')  # probably best to reorganize the repository
    parser.add_argument('--solvent', default='water', help='Name of solvent')
    parser.add_argument('--tau_t', default=0.1, type=float, help='Temperature coupling time constant')
    parser.add_argument('--tau_p', default=20, type=float, help='Pressure coupling time constant')

    args = parser.parse_args()

    return args


class SimulationMdp(object):

    def __init__(self, gro, title='MD Simulation', T=300, em_steps=5000, ensemble='npt', time_step=0.002, length=1000,
                 frames=500, p_coupling='semiisotropic', barostat='berendsen', genvel='yes', restraints=False,
                 xlink=False, bcc=False, tau_p=20, tau_t=0.1):
        """
        :param gro: (str) coordinate file which will be simulated
        :param title: (str) name of simulation
        :param T: (float) simulation temperature
        :param em_steps: (int) number of steps to take during energy minimization
        :param ensemble: (str) thermodynamic ensemble to simulate system in
        :param time_step: (float) simulation time step (fs)
        :param length: (int) simulation length, picoseconds
        :param frames: (int) frequency (number of steps) of recording positions, velocity, energy etc.
        :param p_coupling: (str) type of pressure coupling (None, semiisotropic, isotropic etc.)
        :param barostat: (str) barostat to use for pressure control
        :param genvel: (str) 'yes' if velocities should be generated for initial configuration. 'No' (or anything else)
        if you want to use velocities already present in coordinate file
        :param restraints: (bool) whether or not the system has been restrained (meaning a special topology file has
        been created
        :param xlink: (bool) whether the system is being run through the crosslinking algorithm
        :param bcc: (bool) if we are simulating the bicontinous cubic system
        :param tau_p: (int) time constant for pressure coupling
        :param tau_t: (int) time constant for temperature coupling
        """

        # initialize variables
        self.top = SystemTopology(gro, restraints=restraints)
        self.gro = gro
        self.title = title
        self.temperature = float(T)
        self.em_steps = int(em_steps)
        self.ensemble = ensemble
        self.time_step = float(time_step)
        self.length = int(length)
        self.frames = int(frames)
        self.p_coupling = p_coupling
        self.barostat = barostat
        self.genvel = genvel
        self.restraints = restraints
        self.xlink = xlink
        self.bcc = bcc
        self.tau_p = tau_p
        self.tau_t = tau_t

    def write_em_mdp(self, out='em.mdp'):
        """
        :param out: (str) name of output file
        """

        title = 'title = Energy Minimization'
        integrator = 'integrator = steep'
        nsteps = 'nsteps = %s' % self.em_steps
        cutoff_scheme = 'cutoff-scheme = verlet'
        nstlist = 'nstlist = 40'

        f = open('%s' % out, 'w')
        f.writelines([title + '\n', integrator + '\n', nsteps + '\n', cutoff_scheme + '\n', nstlist + '\n'])

    def write_npt_mdp(self):

        a = []
        a.append(['title = NPT simulation of %s at %s K\n' % (self.gro, self.temperature)])
        # a.append(['cutoff-scheme = verlet'])  # I think verlet is default
        a.append(['integrator = md\n'])  # this also might be default
        a.append(['dt = %s\n' % self.time_step])
        a.append(['nsteps = %s\n' % int(self.length / self.time_step)])
        a.append(['continuation = no\n'])
        a.append(['constraints = h-bonds\n'])
        a.append(['constraint-algorithm = lincs\n'])
        a.append(['cutoff-scheme = Verlet\n'])
        nstx = int(self.length / (self.time_step * self.frames))  # output frequency
        a.append(['nstxout = %s\n' % nstx])
        a.append(['nstvout = %s\n' % nstx])
        a.append(['nstfout = %s\n' % nstx])
        a.append(['nstenergy = %s\n' % nstx])
        a.append(['nstlist = 40\n'])  # potential inputs in the future vvv
        a.append(['nstype = grid\n'])
        a.append(['vdwtype = PME\n'])
        a.append(['coulombtype = PME\n'])
        a.append(['Tcoupl = v-rescale\n'])
        a.append(['tc-grps = system\n'])
        a.append(['tau-t = %s\n' % self.tau_t])
        a.append(['ref-t = %s\n' % self.temperature])
        if not self.restraints:
            a.append(['Pcoupl = %s\n' % self.barostat])
            a.append(['Pcoupltype = %s\n' % self.p_coupling])
            if self.barostat == 'Parrinello-Rahman':
                a.append(['tau-p = 20\n'])  # tau-p should be at least  20 times larger than nstpcouple*dt. nstpcoupl defaults to the value of nstlist (40)
            if self.p_coupling == 'Isotropic':
                a.append(['ref-p = 1\n'])
                a.append(['compressibility = 4.5e-5\n'])
            else:
                a.append(['ref-p = %s\n' % ' '.join([str(1) for i in self.top.residues])])
                a.append(['compressibility = 4.5e-5 4.5e-5\n'])
        if self.genvel == 'yes':
            a.append(['gen-vel = yes\n'])
            a.append(['gen-temp = %s\n' % self.temperature])
        else:
            a.append(['gen-vel = no\n'])
        a.append(['pbc = xyz\n'])
        a.append(['DispCorr = EnerPres\n'])
        if self.xlink:
            a.append('periodic-molecules = yes\n')
            a.append('lincs-iter=4')
        if self.restraints:
            a.append(['refcoord-scaling = all\n'])

        with open('%s.mdp' % self.ensemble, 'w') as f:
            for line in a:
                f.write(line[0])

    def write_nvt_mdp(self):

        a = []
        a.append(['title = NVT simulation of %s\n' % self.gro])
        # a.append(['cutoff-scheme = verlet'])  # I think verlet is default
        a.append(['integrator = md\n'])  # this also might be default
        a.append(['dt = %s\n' % self.time_step])
        a.append(['nsteps = %s\n' % int(self.length / self.time_step)])
        a.append(['continuation = no\n'])
        a.append(['constraints = h-bonds\n'])
        a.append(['constraint-algorithm = lincs\n'])
        nstx = int(self.length / (self.time_step * self.frames))  # output frequency
        a.append(['nstxout = %s\n' % nstx])
        a.append(['nstvout = %s\n' % nstx])
        a.append(['nstfout = %s\n' % nstx])
        a.append(['nstenergy = %s\n' % nstx])
        a.append(['nstlist = 40\n'])
        a.append(['nstype = grid\n'])
        a.append(['vdwtype = PME\n'])
        a.append(['coulombtype = PME\n'])
        a.append(['Tcoupl = v-rescale\n'])
        a.append(['tc_grps = system\n'])
        a.append(['tau_t = %s\n' % self.tau_t])
        a.append(['ref_t = %s\n' % self.temperature])
        a.append(['gen_vel = yes\n'])
        a.append(['pbc = xyz\n'])
        a.append(['DispCorr = EnerPres\n'])
        if self.xlink:
            a.append('periodic-molecules = yes\n')
        if self.restraints:
            a.append(['refcoord_scaling = all\n'])

        with open('%s.mdp' % self.ensemble, 'w') as f:
            for line in a:
                f.write(line[0])


if __name__ == "__main__":

    args = initialize()

    mdp = SimulationMdp(args.gro, title=args.title, T=args.temp, em_steps=args.em_steps, ensemble=args.ensemble,
                        time_step=args.dt, length=args.length, frames=args.frames, p_coupling=args.pcoupltype,
                        barostat=args.barostat, genvel=args.genvel, restraints=args.restraints, xlink=args.xlink,
                        bcc=args.bcc, tau_p=args.tau_p, tau_t=args.tau_t)

    mdp.write_em_mdp()  # write energy minimization .mdp without asking

    if args.ensemble == 'npt':
        mdp.write_npt_mdp()
    elif args.ensemble == 'nvt':
        mdp.write_nvt_mdp()