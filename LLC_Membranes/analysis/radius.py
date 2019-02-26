#!/usr/bin/env python

from LLC_Membranes.analysis.molecular_geometry import Geometry
import sqlite3 as sql
import argparse
import os


def initialize():

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--traj', default='npt_whole.xtc', type=str, help='MD trajectory with whole molecules')
    parser.add_argument('-g', '--gro', default='npt.gro', type=str, help='Gro associated with trajectory')
    parser.add_argument('-r', '--residue', default='ACH', type=str, help='Name of residue whose radius is desired')
    parser.add_argument('-d', '--tablename', default='size', type=str, help='Name of table in database of same name to update')

    return parser


if __name__ == "__main__":

    args = initialize().parse_args() 

    geom = Geometry(args.gro, args.traj, args.residue)
    geom.calculate_radius()
    mean_radius = geom.radius.mean()
    std_radius = geom.radius.std()
    nsolute = geom.nres

    path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    connection = sql.connect('%s/%s.db' % (path, args.tablename))
    crsr = connection.cursor()

    check_existence = "SELECT COUNT(1) FROM %s WHERE name = '%s' and nsolute = %d" %(args.tablename, args.residue, nsolute)

    exists = crsr.execute(check_existence).fetchall()[0][0]

    if exists: 
        fill_new_entry = "UPDATE %s SET end_to_end = %.3f, end_to_end_std = %.3f WHERE name = '%s' and " \
                         "nsolute = %d" % (args.tablename, mean_radius, std_radius, args.residue, nsolute)
        crsr.execute(fill_new_entry)
    else:
        fill_new_entry = "INSERT INTO %s (name, end_to_end, end_to_end_std, nsolute) VALUES ('%s', %.3f, %.3f, %d)" \
                         % (args.tablename, args.residue, mean_radius, std_radius, nsolute)
        crsr.execute(fill_new_entry)

    connection.commit()
    connection.close()

