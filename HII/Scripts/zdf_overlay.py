#! /usr/bin/env python

import zdf
import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import curve_fit


def initialize():

    parser = argparse.ArgumentParser(description='Calculate a correlation function of positions of atoms in the z direction')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Name of trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop calculations')
    parser.add_argument('--skip', default=1, type=int, help='Usage: --skip n . Sample every nth frame')
    parser.add_argument('-a', '--groups', action='append', nargs='+', type=str,
                        help='Name of atoms to calculate z distribution function with respect to. Each group should be'
                             'passed as separate -a flag.')
    parser.add_argument('-bins', default=1000, type=int, help='Number of bins to use when binning distances')
    parser.add_argument('-apl', default=5, type=float, help='Atoms or monomers per layer')
    parser.add_argument('-l', '--load', type=str, help='Load compressed numpy array')
    parser.add_argument('-com', '--com', action="store_true", help='Calculate zdf based on com of atoms')
    parser.add_argument('-r', '--res', type=str, help='Residue name to calculate zdf of')
    parser.add_argument('-cl', type=float, help='Correlation time measured from experiment')

    args = parser.parse_args()

    return args


def grps(atoms):

    groups = []

    for i in atoms:

        if i[0] == 'tails' or i[0] == 'Tails':

            groups.append(['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21',
                 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35',
                 'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48'])

        elif i[0] == 'Head Groups':  # needs to be passed to argparse with quotes around it
            groups.append(['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O', 'O1', 'O2', 'O3', 'O4'])

        else:
            groups.append(i)

    return groups


# def decay(p, x):
#     """
#     :param p: parameters: [period of oscillations, correlation time, amplitude, phase shift]
#     :param x: x axis values
#     :return: exponential function that sinusoidially decaying to one
#     """
#
#     return 1 + p[2]*np.sin(p[0]*x + p[3])*np.exp(-x/p[1])


def decay(x, a, b, c, d):
    """
    :param p: parameters: [period of oscillations, correlation time, amplitude, phase shift]
    :param x: x axis values
    :return: exponential function that sinusoidially decaying to one
    """

    return 1 + c*np.sin(a*x + d)*np.exp(-x/b)


def errorfunc(p, x, z):
    return decay(p, x) - z


if __name__ == "__main__":

    args = initialize()
    groups = grps(args.groups)  # in case of special groups such as 'Tails' or 'Head Groups'

    if args.load:

        zdf_overlay = np.load("%s" % args.load)
        zdf_avg = zdf_overlay["zdf_avg"]
        z = zdf_overlay["z"]
        z_avg = zdf_overlay["z_avg"]

    else:

        groups = grps(args.groups)  # in case of special groups such as 'Tails' or 'Head Groups'

        t = md.load(args.traj, top=args.gro)[args.begin:args.end:args.skip]

        box = t.unitcell_vectors
        L = np.mean(box[:, 2, 2])  # average z length of unit cell

        zdf_avg = np.zeros([len(groups), args.bins])

        for i, grp in enumerate(groups):

            for atom in grp:

                if args.res:
                    keep = [a.index for a in t.topology.atoms if a.name == atom and a.residue.name == args.res]
                else:
                    keep = [a.index for a in t.topology.atoms if a.name == atom]

                pos = t.xyz[:, keep, :]

                periodic = zdf.z_periodic(pos, box)

                z, d = zdf.zdf(periodic, 4, args.apl, box, args.bins, atom)

                zdf_avg[i, :len(d)] += d

            zdf_avg[i, :] /= len(grp)  # average of all zdfs for this group

        z_avg = len(keep)/(4*L)  # divide by 4 since there are 4 pores

        np.savez_compressed("zdf_overlay", zdf_avg=zdf_avg, z=z, z_avg=z_avg)

    plt.figure()

    for i in range(len(groups)):
        plt.plot(z, zdf_avg[i, :len(z)]/z_avg, label=args.groups[i][0])  # divide L by three since there are periodic copies in the +/- x direction

    period = 0.438
    p = np.array([2*np.pi/period, 2, 2, 1])

    start = 0
    while z[start] < 0.5:
        start += 1

    end = -3

    solp, cov_x = curve_fit(decay, z[start:end], zdf_avg[0, start:len(z) + end]/z_avg, p)

    plt.plot(z[start:], decay(z[start:], solp[0], solp[1], solp[2], solp[3]), '--', color='black', label='Least squares fit')

    print('Correlation length = %1.2f +/- %1.2f angstroms' % (10*solp[1], 10*np.sqrt(cov_x[1, 1])))
    print('Oscillation Period = %1.3f nm' % (10*(2*np.pi)/solp[0]))

    if args.cl:
        start = int(0.1*len(z))
        plt.plot(z[start:], 1 + np.exp(-z[start:]/args.cl), '--', color='black', label='Correlation function')

    plt.xlabel('Z distance separation (nm)', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.axes().tick_params(labelsize=14)
    plt.tight_layout()
    plt.legend(loc=1, prop={'size': 16})
    plt.ylim(0, 2)
    plt.tight_layout()
    plt.savefig('zdf_overlay.png')
    plt.show()