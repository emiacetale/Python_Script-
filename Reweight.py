#!/usr/bin/env python

import mdtraj as md
import argparse
from collections import namedtuple
import os 
import subprocess as sp
import numpy as np

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

def write_amber_frame(coords, outfile, comment):
    n_atoms = coords.shape[0]
    print >>outfile, comment
    print >>outfile, '{:d}'.format(n_atoms)
    for xyz0, xyz1 in zip(coords[0::2], coords[1::2]):
        x0, y0, z0 = xyz0*10
        x1, y1, z1 = xyz1*10
	line=crd_format_string.format(x0, y0, z0, x1, y1, z1)
        print >>outfile, line

def energy_GB(filename_a):
    sp.call(["cp",filename_a,"amber.crd"])
    sp.call("pmemd -i amber.in -o mdout -p amber.prmtop -c amber.crd -O", shell=True)
    for line in open('mdinfo'):
        if 'EGB' in line:
	    cols = line.split()
	    ans = float(cols[5])
    return ans  

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

def frame_op(trj, atom_info, basename, resfile):   #iterator over the frames 
    for i, frame in enumerate(trj):       #counter i and frame are now synced
        print "Frame", i 
	#We compute Rham Plot value 
	f=md.compute_phi(frame)
	p=md.compute_psi(frame)
        F=float(f[1][0][0])
	P=float(p[1][0][0])
	#We create the gromacs coord file 
	filename_g = '{}_{:04d}.gro'.format(basename, i)
	#filename="gromacs.gro"
	info=' phi={:8.3f} psi={:8.3f}'.format(F, P)
	#info='test'
	with open(filename_g, 'w') as f:
	    write_gro_frame(frame.xyz[0, ...], atom_info, f, info)
	
	#We also create the Amber coord file 
    filename_a = '{}_{:04d}.inpcrd'.format(basename, i)
    with open(filename_a,'w') as f:
	    write_amber_frame(frame.xyz[0, ...], f, info)
	
	#We compute the SEA and GB Energy
	Egb   = energy_GB(filename_a)
	Esea  = energy_SEA(filename_g)

	#And we write the output in a file 
	line='{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:11.6f}'.format(F, P, Esea, Egb, np.exp(-(Esea-Egb)/0.59616123))
	print >>resfile, line

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('trj_filename', help='trajectory file')
    parser.add_argument('pdb_filename', help='pdb file')
    parser.add_argument('basename', help='otuput coord names')
    parser.add_argument('resfile', help='Result P distr file ')
    return parser.parse_args()

def main():
    args = parse_args()
    trj=md.load(args.trj_filename, top=args.pdb_filename)
    atom_info = get_atom_info(trj)
    with open(args.resfile,"w") as resfile: 
         print >>resfile, "phi psi Esea Egb P(SEA)"
         frame_op(trj, atom_info, args.basename, resfile)
    #frame_op(trj, atom_info, args.basename)


if __name__ == '__main__': #Weird Python way to execute main()
    main()

