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

from simtk.openmm.app import AmberPrmtopFile, OBC2, GBn, GBn2, Simulation, PDBFile, StateDataReporter
from simtk.openmm.app import forcefield as ff
from simtk.openmm import LangevinIntegrator, Platform
from simtk.unit import kelvin, picosecond, femtosecond, angstrom, nanometer, kilojoule_per_mole
from simtk.unit import Quantity, kilojoule, mole, gram

kT= 300*1.9872041E-3  #In kcal/mol

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

def frame_op(sim0, sim1, trj_filename):
    E0=[]
    E1=[]
    for frame in md.iterload(trj_filename, chunk=1):
        E0.append(Eframe(sim0,frame))
        E1.append(Eframe(sim1,frame))
    E0=np.array(E0)
    E1=np.array(E1)
    P=-(E1-E0)/kT
    P=np.exp(P)
    W=np.sum(P)
    pDE, pDE_bins = np.histogram((E1-E0),bins=30,density=True) #For quality control 
    return P/W, [pDE, pDE_bins]

def metatrj(P): #Defines how many times a frame should appear in the output trj 
    nf=len(P)             #number of frames
    i=0                   #counter
    abound=np.zeros(nf,dtype=np.int8)  
    while i<nf:
        frame=np.random.randint(low=0,high=nf)
        if P[frame]>np.random.random():
            abound[frame]=abound[frame]+1
            i=i+1
    return abound

def write_files(trj_filename, fram_abound, name_fout, name_phiout, name_psiout):
    i=0
    Fphi=open(name_phiout,'w')
    Fpsi=open(name_psiout,'w')
    frame_list=[]
    with md.formats.PDBTrajectoryFile(name_fout,mode='w') as out:
         for frame in md.iterload(trj_filename, chunk=1):
             nf=0
             psi, phi = dihedral_calc(frame)
             while nf<fram_abound[i]:
                out.write(frame[0].xyz[0, ...],frame.topology)
                Fphi.write("0.0 "+(str(phi[0]).replace('\n',' '))[1:-1]+"\n") #Similar to VMD out
                Fpsi.write((str(psi[0]).replace('\n',' '))[1:-1]+" 0.0 \n")   #Similar to VMD out
                frame_list.append(i)
                nf=nf+1
             i=i+1
    Fphi.close()
    Fpsi.close()
    return frame_list

def Eframe(sim,frame):
    sim.context.setPositions(frame[0].xyz[0, ...])
    state = sim.context.getState(getEnergy=True)
    return ((state.getPotentialEnergy()/kilojoule_per_mole)*0.239005736)

def make_info(frame_list,P,pDE):
    #Per Frame Probability 
    fig_PF=plt.figure()
    ax1 = fig_PF.add_axes([0.25, 0.15, 0.65, 0.8])
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel(r'$ frame $',fontsize=25)
    ax1.set_ylabel(r'$ P(frame) $',fontsize=25)
    p1,=ax1.plot(np.arange(0,len(P),1),P, color='blue', linewidth=2, linestyle="-")
    fig_PF.savefig('Prob_per_frame.png')
    #Aboundancy of original traj frame in the meta-trj
    fig_FA=plt.figure()
    P_R, eR=np.histogram(frame_list,bins=len(P),normed=False)
    ax1 = fig_FA.add_axes([0.25, 0.15, 0.65, 0.8])
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel('frame',fontsize=20)
    ax1.set_ylabel('N of times frame appear in meta traj',fontsize=20)
    p1,=ax1.plot(np.arange(0,len(P),1),P_R, color='blue', linewidth=2, linestyle="-")
    fig_FA.savefig('Abound_F0_in_metatraj.png') 
    print "Used: "+str(len(list(set(frame_list))))+" frames of the original "+str(len(frame_list))+" frames\n"
    np.savetxt("Frame_prob",P) #Writes the frame probability to a file
    #Difference in energy distributiuon (gives an idea of the error of TD perturbation)
    fig_pDE=plt.figure()
    ax1 = fig_pDE.add_axes([0.25, 0.15, 0.65, 0.8])
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel(r'$ \Delta E $',fontsize=25)
    ax1.set_ylabel(r'$ P(\Delta E) $',fontsize=25)
    p1,=ax1.plot(pDE[1][:-1],pDE[0], color='blue', linewidth=2, linestyle="-")
    NP=pDE[0]*np.exp(-pDE[1][:-1]/kT)
    W=np.sum(NP*(pDE[1][1]-pDE[1][0]))
    p2,=ax1.plot(pDE[1][:-1],NP/W, color='red', linewidth=2, linestyle="-")
    exp=ax1.plot(pDE[1][:-1],np.exp(-pDE[1][:-1]/kT), color='black', linewidth=1, linestyle="--")
    ax1.legend((p1,p2),(r'$P_0(\Delta U)$',r'$P_0(    \Delta U) e^{-\beta \Delta U }$'),loc='upper left')
    fig_pDE.savefig('DE_distr.png')
    #List of frames used (some sort of supercompact metatraj)
    Fflist=open('frames_list','w')
    for f in frame_list:
        Fflist.write(str(f)+" ")
    Fflist.close()

def dihedral_calc(trj):
    ind, phi=md.compute_phi(trj,periodic=False,opt=True)
    ind, psi=md.compute_psi(trj,periodic=False,opt=True)
    return (psi*180/np.pi), (phi*180/np.pi)

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
    sim0, sim1=omm_sys(args.amber_top, args.alph_scal,args.beta_scal)  #Pepare sims wih amap
    P, pDE=frame_op(sim0, sim1, args.trj_filename)                     #Calc frames' norm P 
    fram_abound=metatrj(P) 
    print fram_abound
    frame_list=write_files(args.trj_filename, fram_abound, 'meta_traj.pdb', 'phi.txt', 'psi.txt')
    make_info(frame_list,P,pDE)                                        #Create files useful for 
    print "Done :) "

    ##print "E Stored"
    #W=weight(E,trj.topology.n_residues)
    #print W
    ##simV.step(1)

if __name__ == '__main__': #Weird Python way to execute main()
    main()

