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
import matplotlib.pyplot as plt

from simtk.openmm.app import AmberPrmtopFile, OBC2, GBn, GBn2, Simulation, PDBFile, StateDataReporter
from simtk.openmm.app import forcefield as ff
from simtk.openmm import LangevinIntegrator, Platform
from simtk.unit import kelvin, picosecond, femtosecond, angstrom, nanometer, kilojoule_per_mole
from simtk.unit import Quantity, kilojoule, mole, gram

kT= 300*1.9872041E-3  #In kcal/mol

def energy_GB(frame, simV, simGB):
    simV.context.setPositions(frame.xyz[0, ...])
    simGB.context.setPositions(frame.xyz[0, ...])
    stateV = simV.context.getState(getEnergy=True)
    stateGB=simGB.context.getState(getEnergy=True)
    E_V = stateV.getPotentialEnergy()/kilojoule_per_mole
    E_GB=stateGB.getPotentialEnergy()/kilojoule_per_mole 
    return ((E_GB-E_V)*0.239005736)

def omm_sys(topf, a_s, b_s):
    #Read Amber 
    prmtop = AmberPrmtopFile(topf)
    with open(topf) as top_file:
         top = top_file.read()   #Weird thing needed for amap
    #Create 2 systems
    system0 = prmtop.createSystem(nonbondedMethod=ff.CutoffNonPeriodic,
               nonbondedCutoff=100.*nanometer, constraints=ff.AllBonds, implicitSolvent=OBC2)
    system1 = prmtop.createSystem(nonbondedMethod=ff.CutoffNonPeriodic,
               nonbondedCutoff=100.*nanometer, constraints=ff.AllBonds, implicitSolvent=OBC2)
    #Create 2 integrators
    integrator0 = LangevinIntegrator(300*kelvin, 1/picosecond,0.002*picosecond)
    integrator1 = LangevinIntegrator(300*kelvin, 1/picosecond,0.002*picosecond)
    #Add amap to the systems 
    adder = cmap.CMAPAdder(top,1,1)
    adder.add_to_openmm(system0)
    adder = cmap.CMAPAdder(top,a_s,b_s)
    adder.add_to_openmm(system1)
    #Creates the 2 simulations  
    sim0 = Simulation(prmtop.topology, system0, integrator0)
    sim1 = Simulation(prmtop.topology, system1, integrator1)
    return sim0, sim1

def frame_op(sim0, sim1, trj):
    E0=[]
    E1=[]
    for frame in trj:
        E0.append(Eframe(sim0,frame))
        E1.append(Eframe(sim1,frame))
    #print E0
    #print E1
    E0=np.array(E0)
    E1=np.array(E1)
    P=-(E1-E0)/kT
    P=np.exp(P)
    W=np.sum(P)
    return P/W

def metatrj(trj, P, name_out):
    i=0  
    frames_used=[]
    with md.formats.PDBTrajectoryFile(name_out,mode='w') as out:
        while i<len(trj): #We need i since mdtraj cand perform a len() of a trj that is being created
            #print i
            frame=np.random.randint(low=0,high=len(trj))
            if P[frame]>np.random.random():
               out.write(trj[frame].xyz[0, ...],trj.topology)
               frames_used.append(frame)
               i=i+1
    out.close()     
    return frames_used

def Eframe(sim,frame):
    sim.context.setPositions(frame.xyz[0, ...])
    state = sim.context.getState(getEnergy=True)
    return ((state.getPotentialEnergy()/kilojoule_per_mole)*0.239005736)

def make_info(frame_list,P):
    #Per Frame Probability 
    fig_PF=plt.figure()
    ax1 = fig_PF.add_axes([0.25, 0.15, 0.65, 0.8])
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel(r'$ frame $',fontsize=25)
    ax1.set_ylabel(r'$ P(frame) $',fontsize=25)
    p1,=ax1.plot(np.arange(0,len(P),1),P, color='blue', linewidth=2, linestyle="-")
    fig_PF.savefig('Prob_per_frame.pdf')
    #Aboundancy of original traj frame in the meta-trj
    fig_FA=plt.figure()
    P_R, eR=np.histogram(frame_list,bins=len(P),normed=False)
    ax1 = fig_FA.add_axes([0.25, 0.15, 0.65, 0.8])
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel('frame',fontsize=20)
    ax1.set_ylabel('N of times frame appear in meta traj',fontsize=20)
    p1,=ax1.plot(np.arange(0,len(P),1),P_R, color='blue', linewidth=2, linestyle="-")
    fig_FA.savefig('Abound_F0_in_metatraj.pdf') 

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('trj_filename', help='PDB trajectory file')
    #parser.add_argument('pdb_filename', help='pdb file')
    parser.add_argument('amber_top', help='amber topology')
    parser.add_argument('alph_scal', help='scaling of the alpha helix amap correction', type=float)
    parser.add_argument('beta_scal', help='scaling of the beta sheet amap correction', type=float)
    #parser.add_argument('first_frame',  help='1st frame to process', type=int)
    #parser.add_argument('last_frame',  help='last frame to process', type=int)
    return parser.parse_args()

def main():
    args = parse_args()                                                #Get inline input
    trj=md.load_pdb(args.trj_filename)                                 #load traj
    print "trj loaded"
    sim0, sim1=omm_sys(args.amber_top, args.alph_scal,args.beta_scal)  #Pepare sims wih amap
    P=frame_op(sim0, sim1, trj)                                        #Calc frames' norm P 
    print "Prob calculated"
    frame_list=metatrj(trj, P, 'meta_traj.pdb')                        #Create the traj from prob
    print "Trj created"
    make_info(frame_list,P)                                            #Create files useful for 


    ##print "E Stored"
    #W=weight(E,trj.topology.n_residues)
    #print W
    ##simV.step(1)

if __name__ == '__main__': #Weird Python way to execute main()
    main()
