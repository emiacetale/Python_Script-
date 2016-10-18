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


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-trj',    help='trajectory file', nargs='+') 
    parser.add_argument('-prmtop', help='amber top file' , nargs='+') 
    parser.add_argument('-e1', help='1st atom index for end-to-end dist (AMBER/PDB number!) [2]', type=int, default='2')
    parser.add_argument('-e2', help='2nd atom index for end-to-end dist (AMBER/PDB number!) [109]', type=int, default='109')
    return parser.parse_args()
 
def combine_trj_without_water(trj_names,top_names):
    print "loaded TRJ: ",trj_names
    print "loaded TOP: ",top_names
    trj=md.load_mdcrd(trj_names[0], top_names[0])                          #We load the first traj
    atom_index=[]                                                          #atoms # that will be kept in the new traj
    for el in trj.top.atoms:                                               #We select the atom that are not water
        if str(el.residue)[0:3]!='HOH': atom_index.append(el.index) 
    trj=trj.atom_slice(atom_index)                                         #We keep only the non water atom
    for trj_temp_name, top_temp_name in zip(trj_names[1:],top_names[1:]):  #We concatenate the trj:
        trj_temp=md.load_mdcrd(trj_temp_name, top_temp_name)               #   we load them
        trj_temp=trj_temp.atom_slice(atom_index)                           #   we remove water
        trj.join(trj_temp)                                                 #   we join them to the main trj
    return(trj)


def get_end_to_end(trj, a1, a2):                         #We calc end to end dist
    index=np.array([[a1-1,a2-1]])
    return(md.compute_distances(trj,index))

def get_rg(trj):                                         #func isn't really needed at the moment, but maybe in future(?)  
    return(md.compute_rg(trj))

def compute_phi_psi(trj):                                #computes phi and psi 
    phi=md.compute_phi(trj)
    psi=md.compute_psi(trj)
    return phi[1], psi[1]

def p_norm_2D(values):                                   #norm distr (since norm=true is deprecated)
    P_val, edges = np.histogram(values,bins=100,range=(values.min(),values.max()))
    return (P_val/np.trapz(P_val,x=edges[:-1])), edges[:-1]

def print_on_file_2D(X,Y,Fname,string):                  #print on a file, given an array with X and one with Y
    with open(Fname,'w') as F:
       for x, y in zip(X,Y):
           print >>F,string.format(x,y)

def print_value_raw(value,Fname,formatting):            #prints a matrix of values (lenght of lines can be variable)
    formatting=formatting+"   "
    with open(Fname,'w') as F:
        for el in value: 
            for n in el:
                F.write(formatting.format(n))
            F.write("\n")


def print_out(eed,rg,phi,psi):
    P_eed,e_eed=p_norm_2D(eed)
    print_value_raw(np.transpose(np.append([e_eed],[P_eed],axis=0)),'P_R','{:8.6f}')
    P_rg, e_rg=p_norm_2D(rg)
    print_value_raw(np.transpose(np.append([e_rg],[P_rg],axis=0)),'P_Rg','{:8.6f}')
    print_value_raw(phi,'phi','{:8.6f}')
    print_value_raw(psi,'psi','{:8.6f}')

def main():
    args = parse_args()
    trj=combine_trj_without_water(args.trj,args.prmtop)
    eed=get_end_to_end(trj, args.e1, args.e2)
    rg=get_rg(trj)
    phi, psi = compute_phi_psi(trj)
    #print len(eed), len(rg), len(phi), len(psi)
    print_out(eed,rg,phi,psi)

 
if __name__ == '__main__': #Weird Python way to execute main()
    main()


