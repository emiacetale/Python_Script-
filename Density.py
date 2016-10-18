#!/usr/bin/env python

import numpy as np
import mdtraj as md
import argparse


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-xtc', help='xtc trajectory file')
    parser.add_argument('-pdb', help='pdb file (for topology)')
    parser.add_argument('-np', help='number of grid points in each direction [100]', default=35,type=int)
    parser.add_argument('-out_name', help='name of the output file', default='density.cub')
    parser.add_argument('-density', help='bulk number density of the water model you are using [molecules/nm^3]', default=32.85, type=float)
    return parser.parse_args()

def find_edges(frame, npoints):
    bins=np.arange(0.0,frame.unitcell_lengths[0][0]+frame.unitcell_lengths[0][0]/npoints,frame.unitcell_lengths[0][0]/npoints)
    return (bins,bins,bins)

def get_frame_density(coord,edges):
    d,e=np.histogramdd(coord,bins=edges,normed=False,)
    return d

def print_density(d,edges,file_name,np):
    with open(file_name,'w') as outfile:
         step=edges[0][1]-edges[0][0]
         offset=step/2.0
         print >>outfile,'Comment line 1'
         print >>outfile,'Comment line 2'
         line='{:3d}\t{:12.6f}\t{:12.6f}\t{:12.6f}'.format(0,offset*10*1.88972,offset*10*1.88972,offset*10*1.88972)
         print >>outfile, line
         line='{:3d}\t{:12.6f}\t{:12.6f}\t{:12.6f}'.format(np,step*10*1.88972,0.0,0.0)
         print >>outfile, line
         line='{:3d}\t{:12.6f}\t{:12.6f}\t{:12.6f}'.format(np,0.0,step*10*1.88972,0.0)
         print >>outfile, line
         line='{:3d}\t{:12.6f}\t{:12.6f}\t{:12.6f}'.format(np,0.0,0.0,step*10*1.88972)
         print >>outfile, line
         for plane in d: 
             for line in plane:
                 for point in line:
                     outfile.write('{:12.6f} '.format(point))
                 outfile.write('\n')
             outfile.write('\n')
         

def main():
    args = parse_args()
    #We need to prepare some things
    single_frame=md.load(args.pdb)                            
    edges=find_edges(single_frame, args.np)
    mask=single_frame.topology.select("water and name O")                              #Index of the atom to keep
    density=np.zeros((args.np, args.np,args.np),dtype=np.float)                        #To store results
    nf=0
    for frame in md.iterload(args.xtc, top=args.pdb, chunk=1):                         #Loop over the frames (one by one for mem)
        density += get_frame_density(np.take(frame[0].xyz[0],mask,axis=0),edges)
        nf += 1                                                                        #Any more elegant way for this?  
        #print density 
    density=density/nf
    density=density/(args.density*((edges[0][1]-edges[0][0])*(edges[0][1]-edges[0][0])*(edges[0][1]-edges[0][0])))
    print_density(density,edges,args.out_name,args.np)

if __name__ == '__main__': #Weird Python way to execute main()
    main()

