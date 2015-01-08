#!/usr/bin/env python
 
from simtk.openmm.app import AmberPrmtopFile, OBC2, GBn, GBn2, Simulation,PDBFile, StateDataReporter, PDBReporter
from simtk.openmm.app import forcefield as ff
from simtk.openmm import LangevinIntegrator, MeldForce, Platform, RdcForce, CustomExternalForce
from simtk.unit import kelvin, picosecond, femtosecond, angstrom, nanometer
from simtk.unit import Quantity, kilojoule, mole, gram
from meld.system.restraints import SelectableRestraint, NonSelectableRestraint, DistanceRestraint, TorsionRestraint
from meld.system.openmm_runner import cmap
from meld.system.restraints import ConfinementRestraint, DistProfileRestraint, TorsProfileRestraint, CartesianRestraint
from meld.system.restraints import RdcRestraint

import mdtraj as md
import argparse
from collections import namedtuple
import os 
import subprocess as sp
import numpy as np

import matplotlib as mpl    #To create graph directly in the queue 
mpl.use('Agg')              #"   "     "     "        "  "   " 

import matplotlib.pyplot as plt

from simtk.openmm import app
from simtk.openmm.app import AmberPrmtopFile, OBC2, GBn, GBn2, Simulation, PDBFile, StateDataReporter
from simtk.openmm.app import forcefield as ff
from simtk.openmm import LangevinIntegrator, Platform
from simtk.unit import kelvin, picosecond, femtosecond, angstrom, nanometer, kilojoule_per_mole
from simtk.unit import Quantity, kilojoule, mole, gram
import time

kT= 300*1.9872041E-3  #In kcal/mol

def omm_sys(crdf,topf):          #Prepare the sys
    #Read Amber
    crd = app.AmberInpcrdFile(crdf)
    prmtop = AmberPrmtopFile(topf)
    with open(topf) as top_file:
         top = top_file.read()   #Weird thing needed for amap
    system = prmtop.createSystem(nonbondedMethod=ff.CutoffNonPeriodic, 
           nonbondedCutoff=1.8*nanometer, constraints=ff.AllBonds, implicitSolvent=OBC2) 
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond,0.002*picosecond)
    adder = cmap.CMAPAdder(top,1,1)
    adder.add_to_openmm(system)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    #sim = Simulation(prmtop.topology, system, integrator)
    sim = Simulation(prmtop.topology, system, integrator, platform, properties)
    sim.context.setPositions(crd.positions)
    sim.minimizeEnergy()         #We prepare the sim
    sim.context.setVelocitiesToTemperature(300*kelvin)
    sim.step(100)
    return sim

def benchmarkings(sim, b_time, rep):
    T = []
    C = []
    for n in range(rep):
        t,c=benchmarking(sim, b_time)
        T.append(t)
        C.append(c)
    return np.average(np.array(T)), np.average(np.array(C))

def benchmarking(sim, b_time):
    t0=time.time()
    c0=time.clock()
    sim.step(int(b_time/0.002))
    t1=time.time()
    c1=time.clock()
    return t1-t0, c1-c0

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Prot name (all files should have the same name and different extensions)')
    parser.add_argument('-b_time', help='lenght of Benchmarking time (ps) [50]', default=50.0, type=float)
    parser.add_argument('-n_rep', help='number of time to repeat the sim (for averaging) [5]', default=5, type=int)
    return parser.parse_args()

def main():
    args = parse_args()                            #Get inline input
    simulation = omm_sys(args.filename+'.inpcrd', args.filename+'.prmtop') #Prep the sys
    print 'ready'
    time, clock = benchmarkings(simulation,args.b_time,args.n_rep)      #Perofrm a benchmark 
    print args.filename,': TIME:  ',time
    print args.filename,': CLOCK: ',clock

if __name__ == '__main__': #Weird Python way to execute main()
    main()

