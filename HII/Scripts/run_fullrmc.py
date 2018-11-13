#!/usr/bin/env python2

import os
import numpy as np
import argparse
from fullrmc.Globals import LOGGER, FLOAT_TYPE
from fullrmc.Engine import Engine
from fullrmc.Core.Collection import rebin
from fullrmc.Constraints.DistanceConstraints import InterMolecularDistanceConstraint
from fullrmc.Constraints.StructureFactorConstraints import StructureFactorConstraint
from fullrmc.Core.MoveGenerator import MoveGeneratorCombinator
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Generators.Rotations import RotationGenerator


def initialize():

    parser = argparse.ArgumentParser(description='Add missing fields to .pdb for fullrmc run')

    parser.add_argument('-i', '--pdb', default='initial_elements.pdb', help='Name of .pdb file to modify')
    parser.add_argument('-e', '--exp', default='experimental_normalized.fq', help='Name of file with experimental data')
    parser.add_argument('-o', '--out', default='final.pdb', help='Name of output structure at end of run')
    parser.add_argument('-c', '--cont', action="store_true", help='Use if engine is already made')
    parser.add_argument('-n', '--nsteps', default=100000, type=int, help='Number of monte carlo moves to make')
    parser.add_argument('--cores', default=4, type=int, help='Number of cores to run on')
    parser.add_argument('-ta', '--translation_amplitude', type=float, help='Max translation distance a molecule can jump'
                                                                         '(angstroms)')
    parser.add_argument('-ra', '--rotation_amplitude', type=float, help='Max rotation angle a molecule can turn (degrees)')

    args = parser.parse_args()

    return args


def run_normal(nsteps, saveFrequency, engineFilePath):

    ENGINE.set_groups_as_molecules()
    ENGINE.run(numberOfSteps=nsteps, saveFrequency=saveFrequency, restartPdb=outputpdbFilePath, xyzFrequency=1000, xyzPath=xyzFilePath, ncores=args.cores)

args = initialize()
# Create engine
# DIR_PATH = os.path.dirname( os.path.realpath(__file__) )
DIR_PATH = os.getcwd()

engineFileName = "engine.rmc"
sqFileName = args.exp
pdbFileName = args.pdb
outputpdbFileName = args.out
xyzFileName = "checkpoint.xyz"

# engine variables
sqExpPath = os.path.join(DIR_PATH, sqFileName)
pdbPath = os.path.join(DIR_PATH, pdbFileName)
engineFilePath = os.path.join(DIR_PATH, engineFileName)
outputpdbFilePath = os.path.join(DIR_PATH, outputpdbFileName)
xyzFilePath = os.path.join(DIR_PATH, xyzFileName)

if args.cont:
    FRESH_START = False
else:
    FRESH_START = True

# check if Engine exists. If not build it, otherwise load it.
ENGINE = Engine(path=None)
if not ENGINE.is_engine(engineFilePath) or FRESH_START:
    # create engine
    ENGINE = Engine(path=engineFilePath, freshStart=True)
    ENGINE.set_pdb(str(pdbPath))
    # add S(Q) experimental data
    Sq = np.transpose(rebin(np.loadtxt(sqExpPath), bin = 0.05)).astype(FLOAT_TYPE)
    SF_CONSTRAINT = StructureFactorConstraint(experimentalData=Sq, weighting="atomicNumber")
    ENGINE.add_constraints([SF_CONSTRAINT])
    print('Structure factor constraint added')
    EMD_CONSTRAINT = InterMolecularDistanceConstraint()
    ENGINE.add_constraints([EMD_CONSTRAINT])
    print('Intermolecular Distance constraint added')
    ENGINE.save()
    print('Engine Saved')
else:
    ENGINE = ENGINE.load(engineFilePath)
    SF_CONSTRAINT, EMD_CONSTRAINT = ENGINE.constraints

print('Generating Moves')

for g in ENGINE.groups:
    TMG = TranslationGenerator(amplitude=args.translation_amplitude)
    if len(g) > 1:
        RMG = RotationGenerator(amplitude=args.rotation_amplitude)
        MG = MoveGeneratorCombinator(combination=[TMG, RMG], shuffle=True)
    else:
        MG = MoveGeneratorCombinator(combination=[TMG], shuffle=True)
    g.set_move_generator(MG)

print('Starting Run')
run_normal(nsteps=args.nsteps, saveFrequency=1000, engineFilePath=engineFilePath)
