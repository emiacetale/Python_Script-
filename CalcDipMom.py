#!/usr/bin/env python

import numpy as np
#import matplotlib as mpl    #To create graph directly in the queue 
#mpl.use('Agg')              #"   "     "     "        "  "   " 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time
import matplotlib.font_manager as fm
import mdtraj as md
import argparse

def get_mus(trj,charges):
    mu=[]
    for frame in trj: 
        mu.append(get_mu(frame.xyz,charges))
    return np.array(mu)

def get_mu(coord,q):
    cog=np.average(coord,axis=1)        #the center of geometry is the average of the coord
    #cog=np.array([0.0, 0.0, 0.0])
    r=(coord-cog)                       #we calc the vectors between the cog and every atom
    r=r*q[np.newaxis, :].T             #we multiply for the q "weight" (weird way of numpy to transpose a vector)
    mu_vect=np.sum(r,axis=1)            #the mu vector
    mu=np.linalg.norm(mu_vect)          #the magnitude of the dipole moment
    return mu

def get_q(top_file):
    with open(top_file,'r') as top:
         content=top.read()
    content=content.split('[')
    q=[]
    for sector in content:
       if 'atoms ]' in sector:
           lines=sector.split('\n')
           for line in lines[1:]: 
               line=line.split()
               if len(line)==8: 
                   q.append(float(line[6]))
    return np.array(q)         

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-trj', help='trajectory file')
    parser.add_argument('-pdb', help='pdb file')
    parser.add_argument('-top', help='top file')
    return parser.parse_args()

def main():
    args = parse_args()
    trj=md.load(args.trj, top=args.pdb)
    q=get_q(args.top)
    mu=get_mus(trj,q)
    print mu.mean()
 
if __name__ == '__main__': #Weird Python way to execute main()
    main()

'''
nbins=181
P_dih=np.array([[5,6,19,20]])-1 #The numbers are from the PDB so we subtract 1 to go to pytnon numbering
V_dih=np.array([[5,6,21,22]])-1 #The numbers are from the PDB so we subtract 1 to go to pytnon numbering

traj_P=t = md.load('./PRO/sim/trj_onlyAA.xtc',top='./PRO/PRO.pdb')
traj_V=t = md.load('./VAL/sim/trj_onlyAA.xtc',top='./VAL/VAL.pdb')

theta_P=md.compute_dihedrals(traj_P,P_dih)
theta_V=md.compute_dihedrals(traj_V,V_dih)

dist_P, edges=np.histogram(theta_P,bins=nbins,range=(-np.pi,+np.pi))
dist_V, edges=np.histogram(theta_V,bins=nbins,range=(-np.pi,+np.pi))

dist_P=dist_P/np.trapz(dist_P,x=edges[:-1])   #Norm the distr
dist_V=dist_V/np.trapz(dist_V,x=edges[:-1])   #Norm the distr
#print len(edges), len(dist_P)
fig=plt.figure()
ax1 = fig.add_axes([0.15, 0.2, 0.8, 0.7]) 
p1,=ax1.plot(edges[:-1],dist_P, color='blue', linewidth=2, linestyle="-")
p2,=ax1.plot(edges[:-1],dist_V, color='red', linewidth=2, linestyle="-")

ax1.set_xlim(-np.pi, +np.pi)
ax1.set_ylim(0, 2.0)
l1 = ax1.legend([p1, p2], ["PRO", "VAL"], loc=2, ncol=1, frameon=False, prop=fm.FontProperties(size=25))
ax1.set_xlabel(r'$\theta_{C=O, C=O}$ / rad',fontsize=25)
ax1.set_ylabel(r'$P(\theta) / 10^{-2}$',fontsize=25)
ax1.set_xticks(np.arange(-np.pi,np.pi+0.1,np.pi/2))
ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi$',r'$\pi/2$'])
ax1.set_yticks((0.0, 0.5, 1.0, 1.5, 2.0))
ax1.tick_params(axis='both', which='major', labelsize=20,direction='in', pad=10)

plt.show()
fig.savefig('Dihedral_comparison.png')

'''

