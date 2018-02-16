#!/usr/bin/env python

import numpy as np
import argparse
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Measure Membrane Thickness')

    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help='Path to input file')
    parser.add_argument('-n', '--nframes', default=2, type=int, help='number of frames to make')
    parser.add_argument('-o', '--output', type=str, default='frames', help='Output trajectory')

    args = parser.parse_args()

    return args