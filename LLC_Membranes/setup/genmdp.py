#!/usr/bin/env python

import argparse
from LLC_Membranes.setup.gentop import SystemTopology


def initialize():

    parser = argparse.ArgumentParser(description='Generate .mdp file for simulation of choice')

    parser.add_argument('-T', '--title', default='MD Simulation', type=str, help='Simulation Title')
    parser.add_argument('-g', '--gro', default='initial.gro', type=str, help='coordinate file of system to be simulated')
    parser.add_argument('-t', '--itp', default='dipole.itp', type=str, help='Name of .itp describing monomers')
    parser.add_argument('-s', '--em_steps', default=-1, type=int, help='Steps to take during energy minimization.'
                                                                       'Default is to go forever until convergence.')
    parser.add_argument('-e', '--ensemble', default='npt', type=str, help='Thermodynamic ensemble to put system in')
    parser.add_argument('-d', '--dt', default=0.002, type=float, help='time step (ps)')
    parser.add_argument('-l', '--length', default=1000, type=float, help='simulation length (ps)')
    parser.add_argument('-f', '--frames', default=50, type=int, help='number of frames')
    parser.add_argument('-p', '--pcoupltype', default='semiisotropic', type=str, help='Pressure Couple Type')
    parser.add_argument('--restraints', help='If restraints are on, another mdp option needs to be turned on, so specify '
                                             'this flag', action="store_true")
    parser.add_argument('-x', '--xlink', help='Turn this to "on" if the the system is crosslinked', action="store_true")
    parser.add_argument('-S', '--solvate', help='Specify this if the system has water so an extra line can be added to the '
                                                'topology', action="store_true")
    parser.add_argument('--temp', default=300, help='Specify temperature at which to run simulation')
    parser.add_argument('--mdp', action="store_true", help='Only the .mdp will be written if this option is specified')
    parser.add_argument('--barostat', default='berendsen', type=str, help='pressure coupling scheme to use')
    parser.add_argument('--genvel', default=True, help='generate velocities according to a maxwell'
                                                                     'distribution')
    parser.add_argument('--bcc', action="store_true", help='Generate input files using bicontinuous cubic files')  # probably best to reorganize the repository
    parser.add_argument('--solvent', default='water', help='Name of solvent')
    parser.add_argument('--tau_t', default=0.1, type=float, help='Temperature coupling time constant')
    parser.add_argument('--tau_p', default=20, type=float, help='Pressure coupling time constant')
    parser.add_argument('-nx', '--nstxout', type=int, help='Frequency to output coordinates to trajectory file')
    parser.add_argument('-nv', '--nstvout', type=int, help='Frequency to output velocities to trajectory file')
    parser.add_argument('-nf', '--nstfout', type=int, help='Frequency to output forces to trajectory file')
    parser.add_argument('-ne', '--nstenergy', type=int, help='Frequency to output energy to energy file')

    return parser


class SimulationMdp(object):

    def __init__(self, gro, title='MD Simulation', T=300, em_steps=-1, time_step=0.002, length=1000,
                 p_coupling='semiisotropic', barostat='berendsen', genvel='yes', restraints=False,
                 xlink=False, bcc=False, tau_p=20, tau_t=0.1, nstxout=5000, nstvout=5000, nstfout=5000, nstenergy=5000,
                 frames=None):
        """
        :param gro: (str) coordinate file which will be simulated
        :param title: (str) name of simulation
        :param T: (float) simulation temperature
        :param em_steps: (int) number of steps to take during energy minimization
        :param ensemble: (str) thermodynamic ensemble to simulate system in
        :param time_step: (float) simulation time step (fs)
        :param length: (int) simulation length, picoseconds
        :param p_coupling: (str) type of pressure coupling (None, semiisotropic, isotropic etc.)
        :param barostat: (str) barostat to use for pressure control
        :param genvel: (bool) True if velocities should be generated for initial configuration. \
        if you want to use velocities already present in coordinate file
        :param restraints: (bool) whether or not the system has been restrained (meaning a special topology file has \
        been created
        :param xlink: (bool) whether the system is being run through the crosslinking algorithm
        :param bcc: (bool) if we are simulating the bicontinous cubic system
        :param tau_p: (int) time constant for pressure coupling
        :param tau_t: (int) time constant for temperature coupling
        :param nstxout: frequency of outputting coordinates to trajectory file
        :param nstvout: frequency of outputting velocity to trajectory file
        :param nstfout: frequency of outputting forces to trajectory file
        :param nstenergy: frequency of outputting energy to energy file
        :param frames: number of frames to output. If not None, nstxout, nstvout, nstfou and nstenergy will be \
        adjusted accordingly
        """

        # initialize variables
        self.top = SystemTopology(gro, restraints=restraints, xlink=xlink)
        self.gro = gro
        self.title = title
        self.temperature = float(T)
        self.em_steps = int(em_steps)
        self.time_step = float(time_step)
        self.length = float(length)
        self.p_coupling = p_coupling
        self.barostat = barostat
        self.genvel = genvel
        self.restraints = restraints
        self.xlink = xlink
        self.bcc = bcc
        self.tau_p = tau_p
        self.tau_t = tau_t
        self.em_mdp_name = None
        self.npt_mdp_name = None
        self.nvt_mdp_name = None
        self.nve_mdp_name = None
        self.nstxout = nstxout
        self.nstvout = nstvout
        self.nstfout = nstfout
        self.nstenergy = nstenergy

        if frames is not None:
            self.nstxout = int(self.length / (self.time_step * frames))
            self.nstvout = int(self.length / (self.time_step * frames))
            self.nstfout = int(self.length / (self.time_step * frames))
            self.nstenergy = int(self.length / (self.time_step * frames))

    def write_em_mdp(self, out='em'):
        """
        :param out: (str) name of output file
        """

        title = 'title = Energy Minimization'
        integrator = 'integrator = steep'
        nsteps = 'nsteps = %s' % self.em_steps
        cutoff_scheme = 'cutoff-scheme = verlet'
        nstlist = 'nstlist = 40'
        constraints = 'constraints = h-bonds'

        f = open('%s.mdp' % out, 'w')
        f.writelines([title + '\n', integrator + '\n', nsteps + '\n', cutoff_scheme + '\n', nstlist + '\n', constraints
                      + '\n'])

        if self.xlink:
            f.write('periodic-molecules = yes\n')

        self.em_mdp_name = out

    def write_npt_mdp(self, out='npt'):
        """
        :param out: (str) name of output file
        """

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
        a.append(['nstxout = %s\n' % self.nstxout])
        a.append(['nstvout = %s\n' % self.nstvout])
        a.append(['nstfout = %s\n' % self.nstfout])
        a.append(['nstenergy = %s\n' % self.nstenergy])
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
                a.append(['ref-p = %s\n' % ' '.join(['1', '1'])])
                a.append(['compressibility = 4.5e-5 4.5e-5\n'])
        if self.genvel:
            a.append(['gen-vel = yes\n'])
            a.append(['gen-temp = %s\n' % self.temperature])
        else:
            a.append(['gen-vel = no\n'])
        a.append(['pbc = xyz\n'])
        a.append(['DispCorr = EnerPres\n'])
        if self.xlink:
            a.append(['periodic-molecules = yes\n'])
            a.append(['lincs-iter=4'])  # I don't think this is necessary
        if self.restraints:
            a.append(['refcoord-scaling = all\n'])

        with open('%s.mdp' % out, 'w') as f:
            for line in a:
                f.write(line[0])

        self.npt_mdp_name = "%s.mdp" % out

    def write_nvt_mdp(self, out='nvt'):
        """
        :param out: (str) name of output file
        """
        a = []
        a.append(['title = NVT simulation of %s\n' % self.gro])
        # a.append(['cutoff-scheme = verlet'])  # I think verlet is default
        a.append(['integrator = md\n'])  # this also might be default
        a.append(['dt = %s\n' % self.time_step])
        a.append(['nsteps = %s\n' % int(self.length / self.time_step)])
        a.append(['continuation = no\n'])
        a.append(['constraints = h-bonds\n'])
        a.append(['constraint-algorithm = lincs\n'])
        a.append(['nstxout = %s\n' % self.nstxout])
        a.append(['nstvout = %s\n' % self.nstvout])
        a.append(['nstfout = %s\n' % self.nstfout])
        a.append(['nstenergy = %s\n' % self.nstenergy])
        a.append(['nstlist = 40\n'])
        a.append(['nstype = grid\n'])
        a.append(['vdwtype = PME\n'])
        a.append(['coulombtype = PME\n'])
        a.append(['Tcoupl = v-rescale\n'])
        a.append(['tc_grps = system\n'])
        a.append(['tau_t = %s\n' % self.tau_t])
        a.append(['ref_t = %s\n' % self.temperature])
        if self.genvel:
            a.append(['gen-vel = yes\n'])
            a.append(['gen-temp = %s\n' % self.temperature])
        else:
            a.append(['gen-vel = no\n'])
        a.append(['pbc = xyz\n'])
        a.append(['DispCorr = EnerPres\n'])
        if self.xlink:
            a.append(['periodic-molecules = yes\n'])
        if self.restraints:
            a.append(['refcoord_scaling = all\n'])

        with open('%s.mdp' % out, 'w') as f:
            for line in a:
                f.write(line[0])

        self.nvt_mdp_name = "%s.mdp" % out

    def write_nve_mdp(self, out='nve'):
        """
        :param out: (str) name of output file
        """
        a = []
        a.append(['title = NVE simulation of %s\n' % self.gro])
        # a.append(['cutoff-scheme = verlet'])  # I think verlet is default
        a.append(['integrator = md\n'])  # this also might be default
        a.append(['dt = %s\n' % self.time_step])
        a.append(['nsteps = %s\n' % int(self.length / self.time_step)])
        a.append(['continuation = no\n'])
        a.append(['constraints = h-bonds\n'])
        a.append(['constraint-algorithm = lincs\n'])
        a.append(['nstxout = %s\n' % self.nstxout])
        a.append(['nstvout = %s\n' % self.nstvout])
        a.append(['nstfout = %s\n' % self.nstfout])
        a.append(['nstenergy = %s\n' % self.nstenergy])
        a.append(['nstlist = 40\n'])
        a.append(['nstype = grid\n'])
        a.append(['vdwtype = PME\n'])
        a.append(['coulombtype = PME\n'])
        if self.genvel:
            a.append(['gen-vel = yes\n'])
            a.append(['gen-temp = %s\n' % self.temperature])
        else:
            a.append(['gen-vel = no\n'])
        a.append(['pbc = xyz\n'])
        a.append(['DispCorr = EnerPres\n'])
        if self.xlink:
            a.append('periodic-molecules = yes\n')
        if self.restraints:
            a.append(['refcoord_scaling = all\n'])

        with open('%s.mdp' % out, 'w') as f:
            for line in a:
                f.write(line[0])

        self.nve_mdp_name = "%s.mdp" % out

    def add_pull_groups(self, ref_groups, coord_groups, k, rate, mdp, geometry='distance', type='umbrella', dim='z'):
        """ Add pull groups

        NOTE: This assumes all options apply to all pull coords (for now)

        :param ref_groups: name of groups used as reference com
        :param coord_groups: name of groups which will be used with a pull coordinate
        :param k: force constant (kJ / mol / nm^2)
        :param rate: pull rate (nm/ps)
        :param mdp: name of .mdp file to add pull parameters to
        :param geometry: how to pull (see http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html)
        :param type: see http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html
        :param dim: axis a long which to pull

        :type ref_groups: list or tuple
        :type coord_groups: list or tuple
        :type k: float
        :type rate: float
        :type mdp: str
        :type geometry: str
        :type type: str
        :type dim: str
        """

        n = len(coord_groups)
        pulldim = ['N', 'N', 'N']
        if 'x' in dim:
            pulldim[0] = 'Y'
        if 'y' in dim:
            pulldim[1] = 'Y'
        if 'z' in dim:
            pulldim[2] = 'Y'
        dim = " ".join(pulldim)

        with open(mdp, "a") as f:
            f.write('\n; Pull code\n')
            f.write('pull = yes\n')
            f.write('pull-ngroups = %d\n' % (len(ref_groups) + len(coord_groups)))
            f.write('pull-ncoords = %d\n' % len(coord_groups))
            f.write('pull-print-components = yes\n')  # print distance from com with sign
            for i, x in enumerate(coord_groups):
                num = i + 1
                f.write('pull-group%d-name = %s\n' % (num, x))
                f.write('pull-coord%d-type = %s\n' % (num, type))
                f.write('pull-coord%d-geometry = %s\n' % (num, geometry))
                f.write('pull-coord%d-dim = %s\n' % (num, dim))
                f.write('pull-coord%d-groups = %d %d\n' % (num, n + i + 1, num))
                f.write('pull-coord%d-rate = %.1f\n' % (num, rate))  # pull rate
                f.write('pull-coord%d-k = %.1f\n' % (num, k))  # set harmonic force constant
                f.write('pull-coord%d-start = yes\n' % num)  # important so that com distance is maintained
                f.write('\n')  # space between groups for clarity
            for i, x in enumerate(ref_groups):
                f.write('pull-group%d-name = %s\n' % (n + i + 1, x))


if __name__ == "__main__":

    args = initialize().parse_args()

    # get output frequencies (important for controlling size of trajectory)
    if args.nstxout:
        nstxout = args.nstxout
    else:
        nstxout = int(args.length / (args.dt * args.frames))

    if args.nstvout:
        nstvout = args.nstvout
    else:
        nstvout = int(args.length / (args.dt * args.frames))

    if args.nstfout:
        nstfout = args.nstfout
    else:
        nstfout = int(args.length / (args.dt * args.frames))

    if args.nstenergy:
        nstenergy = args.nstenergy
    else:
        nstenergy = int(args.length / (args.dt * args.frames))  # output frequency

    mdp = SimulationMdp(args.gro, title=args.title, T=args.temp, em_steps=args.em_steps,
                        time_step=args.dt, length=args.length, p_coupling=args.pcoupltype,
                        barostat=args.barostat, genvel=args.genvel, restraints=args.restraints, xlink=args.xlink,
                        bcc=args.bcc, tau_p=args.tau_p, tau_t=args.tau_t, nstxout=nstxout, nstvout=nstvout,
                        nstfout=nstfout, nstenergy=nstenergy)

    mdp.write_em_mdp()  # write energy minimization .mdp without asking

    if args.ensemble == 'npt':
        mdp.write_npt_mdp()
    elif args.ensemble == 'nvt':
        mdp.write_nvt_mdp()
