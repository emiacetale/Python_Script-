#!/usr/bin/env python

import mdtraj as md
import argparse
import numpy as np
from scipy.spatial import Voronoi,Delaunay

#def pint_res(dg_all,qrange,srange, dgt):
#    DG=dg_all.reshape(len(qrange),len(srange)) # { {DG_s0, DG_s1, ... }@q0 , {DG_s0, DG_s1, ... }@q1 , ... }
#    with open("results",'w') as out:
#      print >>out, "#q_scale  s_scale  DG*   DG*-target"
#      for q, el in zip(qrange, DG):
#          for s, dg in zip(srange, el):
#	      print >>out ,'{:4.2f} {:4.2f} {:9.6f} {:9.6f}'.format( q, s, dg, dg-dgt)   

def find_protein_atoms(pdb):
    t=md.load(pdb)
    return(t.top.select('not water'))

def tetravol(a,b,c,d):
    '''Calculates the volume of a tetrahedron, given vertices a,b,c and d (triplets)'''
    tetravol=abs(np.dot((a-d),np.cross((b-d),(c-d))))/6
    return tetravol


def vol(vor,p):
    '''Calculate volume of 3d Voronoi cell based on point p. Voronoi diagram is passed in v.'''
    dpoints=[]
    vol=0
    for v in vor.regions[vor.point_region[p]]:
        dpoints.append(list(vor.vertices[v]))
    tri=Delaunay(np.array(dpoints))
    for simplex in tri.simplices:
        vol+=tetravol(np.array(dpoints[simplex[0]]),np.array(dpoints[simplex[1]]),np.array(dpoints[simplex[2]]),np.array(dpoints[simplex[3]]))
    return vol

def V_frames_in_trj(traj,topology,atom_index):
    V=[]
    for frame in md.iterload(traj, top=topology,chunk=1):                                                                 
        voronoi = Voronoi(frame.xyz[0])
        volume=0.0
        for i in atom_index:
            volume=volume+vol(voronoi,i)
        V.append(volume)
    return np.array(V)

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-traj', type=str, help='Gromacs XTC file', default='trjfull.xtc')
    parser.add_argument('-pdb', type=str, help='pdb file to use as topology', default='conf.pdb')
    return parser.parse_args()

def main():
    #Script:
    args = parse_args() 
    atom_index = find_protein_atoms(args.pdb)
    V=V_frames_in_trj(args.traj,args.pdb,atom_index)
    
    print '<V>         = '+str(np.average(V))    
    st_dev=np.std(V, dtype=np.float64) 
    print '<V^2>-<V>   = '+str(st_dev)
    print '<V^2>-<V>^2 = '+str(st_dev*st_dev)
    
    #V=[]
    #for frame in md.iterload('output.xtc', top='topology.pdb',chunk=1):
    #    frame = md.load_frame(args.traj, 1, top=args.pdb)
    #    voronoi = Voronoi(frame.xyz[0])
    #    volume=0.0
    #    for i in atom_index:
    #        volume=volume+vol(voronoi,i)
    #    volume

if __name__ == '__main__': #Weird Python way to execute main()
    main()



