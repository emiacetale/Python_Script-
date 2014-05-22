#!/usr/bin/env python

import mdtraj as md
import argparse
from collections import namedtuple
import os 
import subprocess as sp
import numpy as np

from simtk.openmm.app import AmberPrmtopFile, OBC2, GBn, GBn2, Simulation, PDBFile, StateDataReporter
from simtk.openmm.app import forcefield as ff
from simtk.openmm import LangevinIntegrator, Platform
from simtk.unit import kelvin, picosecond, femtosecond, angstrom, nanometer, kilojoule_per_mole
from simtk.unit import Quantity, kilojoule, mole, gram

kT= 300*1.9872041E-3  #In kcal/mol

AtomInfo = namedtuple('AtomInfo', 'res_num res_name atom_num atom_name') #Create a named touple 
gro_format_string = '{:5d}{:5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'          #Define Format
crd_format_string = '{:12.7f}{:12.7f}{:12.7f}{:12.7f}{:12.7f}{:12.7f}'

def get_atom_info(trj):                  #Estract the info we need (A & Res name, and A & R index)                          
    atoms = []
    for atom in trj.topology.atoms:
        #print dir(atom)                 #Show the attributes of atom
        res_num = atom.residue.resSeq
        res_name = atom.residue.name
        atom_num = atom.index
        atom_name = atom.name
        atoms.append(AtomInfo(res_num, res_name, atom_num, atom_name))
    return atoms

def write_gro_frame(coords, atom_info, outfile, comment):
    n_atoms = coords.shape[0]
    assert n_atoms == len(atom_info)
    print >>outfile, comment
    print >>outfile, '{:d}'.format(n_atoms)
    for xyz, atom in zip(coords, atom_info):    #2 syncronized iterator for 2 lists
        x, y, z = xyz
	line = gro_format_string.format(atom.res_num, atom.res_name, atom.atom_name, atom.atom_num+1, x, y, z)
	print >>outfile, line
    print >>outfile, ' 10.000   10.000   10.000' #Fake box size

def energy_GB(frame, simV, simGB):
    simV.context.setPositions(frame.xyz[0, ...])
    simGB.context.setPositions(frame.xyz[0, ...])
    stateV = simV.context.getState(getEnergy=True)
    stateGB=simGB.context.getState(getEnergy=True)
    E_V = stateV.getPotentialEnergy()/kilojoule_per_mole
    E_GB=stateGB.getPotentialEnergy()/kilojoule_per_mole 
    return ((E_GB-E_V)*0.239005736)

def energy_SEA(filename_g):
    sp.call(["cp",filename_g,"gromacs.gro"])
    p=sp.Popen("/home/ebrini/software/SEA_LIBO/bin/solvate -s gromacs -ce none 2>/dev/null", shell=True, stdout=sp.PIPE)
    p.communicate()
    p=sp.Popen("/home/ebrini/software/FSEA/FSEA_adp_big_mol.exe < surface.povdat", shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    #print output
    Ep=float(output)
    p=sp.Popen("rm surface.povdat", shell=True)
    p.communicate()
    return (Ep/4.184) #Our UOM is kcal wile Libo's FSEA is in KJ

def make_files(trj, atom_info, basename):
    info=' File created with a script by EMI' 
    gro_f=[]
    crd_f=[]
    for i, frame in enumerate(trj):
        filename_g = '{}_{:06d}.gro'.format(basename, i)
        gro_f.append(filename_g)
        with open(filename_g, 'w') as f:
             write_gro_frame(frame.xyz[0, ...], atom_info, f, info)
    return gro_f

def omm_sys(top, pdb):
    prmtop = AmberPrmtopFile(top)
    pdb = PDBFile(pdb)
    systemOBC = prmtop.createSystem(nonbondedMethod=ff.CutoffNonPeriodic,
        nonbondedCutoff=100.*nanometer, constraints=ff.HBonds, implicitSolvent=OBC2)
    systemVAC=system = prmtop.createSystem(nonbondedMethod=ff.CutoffNonPeriodic,
               nonbondedCutoff=100.*nanometer, constraints=ff.HBonds, implicitSolvent=None)
    integratorGB = LangevinIntegrator(300*kelvin, 1/picosecond,0.002*picosecond)
    integratorV = LangevinIntegrator(300*kelvin, 1/picosecond,0.002*picosecond)
    simGB = Simulation(prmtop.topology, systemOBC, integratorGB)
    simV  = Simulation(prmtop.topology, systemVAC, integratorV)
    simGB.context.setPositions(pdb.positions)
    simV.context.setPositions(pdb.positions)
    return simV, simGB

def frame_op(simV, simGB, trj, gro_f):
    #print dir(stateV)
    E=[]
    for gf, frame in zip(gro_f,trj):
        #print gf 
        E_SEA=energy_SEA(gf)
        E_GB=energy_GB(frame, simV, simGB)
        E.append([E_SEA,E_GB])
    return E

def weight(E,nr):
    E=np.array(E)
    E=E/nr               #We divide by the number of AA 
    de=np.diff(E,axis=1) #Array of Egb-Esea == -(Esea-Egb)
    de=de/kT               
    de=np.exp(de)        #Array of e^(-dE/kT) == weights
    w=np.sum(de)         #Sum of the weights == cluster weight
    return (w)

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('trj_filename', help='trajectory file')
    parser.add_argument('pdb_filename', help='pdb file')
    parser.add_argument('amber_top', help='amber topology')
    parser.add_argument('first_frame',  help='1st frame to process', type=int)
    parser.add_argument('last_frame',  help='last frame to process', type=int)
    return parser.parse_args()

def cut_trj(trj, first, last):
    if len(trj)<first:
        #print "Not enough frame in the trajectory! (first frame to analyze is > than the last trj frame"
        exit()
    if len(trj)<last:
        last=len(trj)
    return(trj[first:last])


def main():
    #Variables that don't require user input 
    basename='struct_temp' 
    resname='Weights'
    #Actual script 
    args = parse_args()                                         #Get inline input
    trj=md.load_mdcrd(args.trj_filename, top=args.pdb_filename) #load traj
    trj = cut_trj(trj, args.first_frame, args.last_frame)
    atom_info = get_atom_info(trj)                              #get info of the atoms 
    gro_f = make_files(trj, atom_info, basename)                #Make gro file from the traj 
    simV , simGB=omm_sys(args.amber_top, args.pdb_filename)     #makes openMM systems
    E=frame_op(simV, simGB, trj, gro_f)
    #print "E Stored"
    W=weight(E,trj.topology.n_residues)
    print W
    #simV.step(1)

if __name__ == '__main__': #Weird Python way to execute main()
    main()

