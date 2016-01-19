#!/usr/bin/env python

import mdtraj as md
import argparse
import os 
import numpy as np
import matplotlib as mpl    #To create graph directly in the queue 
mpl.use('Agg')              #"   "     "     "        "  "   " 
import matplotlib.pyplot as plt  
import cPickle as pickle

def trim_traj(trj, atom_name, b, e):
    #print b, e 
    index=trj.top.select('name '+atom_name+' and resid '+str(b)+' to '+str(e))
    return(trj.atom_slice(index))

def atom_index(trj, atom_name, b, e):
    index=trj.top.select('name '+atom_name+' and resid '+str(b)+' to '+str(e))
    return(index)

def clean_names(names):
    return [el.split('/')[-1] for el in names] 

def plot_results(dat, names, out_name):
    names=clean_names(names)
    lines=[]
    for d in dat:
        tmp,=plt.plot(d,'-')
        lines.append(tmp)
    plt.legend(lines,names)
    plt.savefig(out_name+'.png')
    #plt.show()

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-traj', help='pdb trajectory file',default='XXX')
    parser.add_argument('-traj_dcd', help='topolgy and trajectory',default='XXX')
    parser.add_argument('-top_dcd', help='topolgy and trajectory',default='XXX')
    parser.add_argument('-refs', help='pdb files to compare', nargs='+')
    parser.add_argument('-atom_name', help='atom name to use for the RMSD calc', default='CA')
    parser.add_argument('-start_res', help='ignore residue before this', type=int, default=8)
    parser.add_argument('-end_res', help='ignore residue after this', type=int, default=52)
    parser.add_argument('-name', help='output name', default='RMSD')
    parser.add_argument('-save_pdb', help='save the pdbs of the atoms used for the RMSD', default=False, action='store_true')
    parser.add_argument('-debug', help='print some info along the way', default=False, action='store_true')
    return parser.parse_args()

def main():
    args = parse_args()
    refs=[]
    pdb_index=[]
    for r in args.refs: 
        refs.append(trim_traj(md.load(r), args.atom_name, args.start_res, args.end_res))

    if args.traj!='XXX':
        print 'Reading pdb trajectory: '+args.traj
        #trj=md.load(args.traj)
        #index_trj=atom_index(traj, args.atom_name, args.start_res, args.end_res)
        trj=trim_traj(md.load(args.traj),args.atom_name, args.start_res, args.end_res)
    if args.traj_dcd!='XXX':
        print 'Reading dcd trajectory: '+args.traj_dcd+' , with topolgy: '+args.top_dcd
        #trj=md.load(args.traj_dcd,top=args.top_dcd)
        #index_trj=atom_index(trj, args.atom_name, args.start_res, args.end_res)
        trj=trim_traj(md.load(args.traj_dcd,top=args.top_dcd),args.atom_name, args.start_res, args.end_res)
    rmsd=[]
    for i,el in enumerate(refs):
        print 'RMSD '+str(i)
        rmsd.append(md.rmsd(trj,el[0],parallel=False, precentered=False))
        if args.debug:
           for a,b in zip(trj.top.residues,el.top.residues):
               print a,b
           for a,b in zip(trj.top.atoms,el.top.atoms):
               print a,b

    #for i,el in enumerate(zip(refs,pdb_index)):
    #    print 'RMSD '+str(i)
    #    rmsd.append(md.rmsd(trj,el[0],parallel=False, atom_indices=index_trj, ref_atom_indices=el[1], precentered=False))
    #print 'plotting'
    plot_results(rmsd, args.refs, args.name)
    with open(args.name+'.pckl','wb') as out:
        pickle.dump(rmsd,out)

    if args.save_pdb:
        trj.save_pdb('traj_'+args.name+'.pdb')
        for i,el in enumerate(refs):
           el.save_pdb('traj_ref_'+args.name+'_'+str(i)+'.pdb')


if __name__ == '__main__': #Weird Python way to execute main()
    main()

