#!/usr/bin/python

import MDAnalysis
import numpy as np
import pickle
import sys
import os

class comp_exp(object):
    """
    Class to calculate the Z dimension bounds of the membrane based on the narrowest
    bounds that satisfy the densityThreshold over all sampled frames using Gromacs.

    Arguments:
      *name*
        full name (including directory) of trr and tpr files
        ASSUMPTION: trr and tpr files have the same name
      *t0*
        first frame of trr to analyze (ps)
      *tf*
        last frame of the trr to analyze (ps)
      *bounds*
        tuple (lower,upper) containing z-dimension bounds of the membrane (nm)
      *len_side*
        box dimension in the x and y dimension (nm)
      *sol_buff*
        distance from z-dimension box ends with sol phase (nm)
        IMPLICIT ASSUMPTION: membrane is roughly centered in z-dimension
      *outdir*
        full output directory (including '/' at the end)
      *outname*
        name of output pickle file
      *dt_restart*
        sampling interval (ps)
    """

    def __init__(self,name,t0,tf,bounds,len_side,sol_buff,outdir='./',outname='comp_exp',
                 dt_restart=100):
        self.name = name
        self.t0 = int(t0)
        self.tf = int(tf)
        b1 = float(bounds[0])
        b2 = float(bounds[1])
        self.bounds = (b1,b2)
        self.len_side = float(len_side)
        self.sol_buff = float(sol_buff)
        self.outdir = outdir
        self.outname = outname
        self.dt_restart = int(dt_restart)

    def _init_arrays(self,t0,tf,dt_restart):
        """
        Initializes density, waterUptake, Cmx, and Csx arrays
        """
        array = np.zeros((tf-t0)/dt_restart+1)
        return array, array.copy(), array.copy(), array.copy()

    def _get_sel(self,u,bounds,sol_buff):
        """
        Returns membrane, sol, poreWaters, memCl, and solCl selections
        """
        mem = u.select_atoms("prop z >= "+str(bounds[0]*10)+" and prop z <= "
                             +str(bounds[1]*10))
        sol = u.select_atoms("prop z <= "+str(sol_buff*10)+" or prop z >= "+
                             str(u.trajectory.ts.dimensions[2]-sol_buff*10))
        poreWaters = mem.select_atoms("resname SOL")
        memCl = mem.select_atoms("resname CL")
        solCl = sol.select_atoms("resname CL")
        return mem, sol, poreWaters, memCl, solCl

    def run(self,**kwargs):
        print("Loading Universe...")
        u = MDAnalysis.Universe(self.name+'.tpr',self.name+'.trr')
        print('Complete!\n')

        density, waterUptake, Cmx, Csx = self._init_arrays(self.t0,self.tf,self.dt_restart)
        memVol = self.len_side**2 * (self.bounds[1]-self.bounds[0]) # nm^3
        solVol = self.len_side**2 * 2 * self.sol_buff # nm^3
        DENSITY_CONVERSION = 1.660539 # conversion from amu/nm^3 to kg/m^3
        AVAGADROS = 6.0221409 * 10 **23 # atoms/mol

        counter = 0
        for ts in u.trajectory[self.t0:self.tf+1:self.dt_restart]:
             mem, sol, poreWaters, memCl, solCl = self._get_sel(u,self.bounds,self.sol_buff)

             density[counter]=mem.total_mass()/memVol * DENSITY_CONVERSION
             waterUptake[counter]=poreWaters.total_mass()/(mem.total_mass()-poreWaters.total_mass())
             Cmx[counter]=float(memCl.n_atoms*10**24)/(memVol*AVAGADROS)
             Csx[counter]=float(solCl.n_atoms*10**24)/(solVol*AVAGADROS)
             counter+=1
             if ts.frame % 500 == 0:
                 print("Analyzing frame: "+str(ts.frame)+"\n")
                 print("CL in membrane: "+str(memCl.n_atoms))

        pickle.dump([density,waterUptake,Cmx,Csx],open(self.outdir+self.outname+'.p','wb'))

def safe_outdir(outdir):
    """
    Creates outdir and all necessary intermediate directories if they don't exist
    """
    dirs = outdir.split('/')
    for i in range(len(dirs)):
        dir_str = '/'
    for j in range(i+1):
        dir_str += dirs[j]+'/'
    if not os.path.exists(dir_str):
        os.system('mkdir '+dir_str)
    return

if __name__ == "__main__":
    name = sys.argv[1]
    t0 = sys.argv[2]
    tf = sys.argv[3]
    bounds = tuple(sys.argv[4]) # nm
    print bounds
    len_side = sys.argv[5] # nm
    sol_buff = sys.argv[6] # nm
    # these need to be specified from rootdir and end in '/'
    outdir = sys.argv[7]
    outname = sys.argv[8]

    # run calculation
    safe_outdir(outdir)
    this_comp_exp = comp_exp(name,t0,tf,bounds,len_side,sol_buff,outdir,outname)
    this_comp_exp.run()