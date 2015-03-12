#!/usr/bin/env python

#import mdtraj as md
import argparse
import numpy as np
import subprocess as sp
import matplotlib as mpl    #To create graph directly in the queue 
mpl.use('Agg')              #"   "     "     "        "  "   " 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time 

def ESEA(traj, top, tpr, new_gro_name):
    gro_files=create_gro(new_gro_name, 'non-Water', traj, tpr) #Create gro files
    E=Eone(gro_files, top)
    return (E*4.184) 

def energy_SEA(filename_g):
    sp.call(["mv",filename_g,"gromacs.gro"])
    p=sp.Popen("/home/ebrini/software/SEA/bin/solvate -s gromacs -ce none -d 12 -i 500 2> /dev/null ", shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    output=output.split('\n')
    for line in output:
        if 'Total' in line: 
            line=line.split()
            E=float(line[1])
    return E 

def Eone(gro_files, top_file):
    sp.call(["cp",top_file,"gromacs.top"])
    E=[]
    for gro_f in gro_files:
         E.append(energy_SEA(gro_f))
    return np.array(E)

def EGromacs(traj, top, grompp, tpr,group):
    p=sp.Popen('trjconv -f '+traj+' -o traj0.trr -s '+tpr+' -pbc mol &> /dev/null ', shell=True, stdin=sp.PIPE)
    p.stdin.write(group+' \n') #We should be sure to remove water from the traj
    p.communicate()[0]            
    p.stdin.close()
    p=sp.Popen('grompp -f '+grompp+' -p '+top+' -c gromacs.gro &> /dev/null ', shell=True, stdin=sp.PIPE)
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('mdrun -rerun traj0.trr -nt 1 &> /dev/null ', shell=True, stdin=sp.PIPE)
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('g_energy -xvg none -o energy.xvg &> /dev/null ', shell=True, stdin=sp.PIPE)
    p.stdin.write('Potential \n')
    p.stdin.write('\n')
    p.communicate()[0]
    p.stdin.close()
    X, E=np.loadtxt('energy.xvg', dtype=float, unpack=True)
    p=sp.Popen('rm energy.xvg ener.edr md.log mdout.mdp topol.tpr traj.trr traj0.trr', shell=True, stdin=sp.PIPE)
    p.communicate()[0]
    p.stdin.close()
    return E
    
def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-trj_l', type=str, help='gromacs trajectory of the solute in solution to read')
    parser.add_argument('-trj_v', type=str, help='gromacs trajectory of the solute in vacuum to read')
    parser.add_argument('-top', type=str, help='gromacs topology to read')
    parser.add_argument('-tpr_l', type=str, help='gromacs tpr to read for solution sim')
    parser.add_argument('-tpr_v', type=str, help='gromacs tpr to read for vac sim')
    parser.add_argument('-qs', type=float, help='scaling of the charge with respect to the original FF')
    parser.add_argument('-ss', type=float, help='scaling of the sigma  with respect to the original FF')
    return parser.parse_args()

def main():
    #"Hard-coded" variables: not much sense to change them at this point, but it is nice to have an handle for future 
    new_top_name="TOP.top"     # Name of topology file 
    new_gro_name="struct"  # Name of conf file 
    gompp_name="grompp.mdp"

    #Script:
    args = parse_args()                                                                      #We read some input 
    Th, Tc, Tt, At, A=read_top(args.top)                                                        #Read and decompose the 
    ATinfo, sigma, epsilon  = get_s_e(At)                                                    #    topology
    Ainfo, q=get_q(A)                                                                        #    ...
    top_file=write_top(Th, Tc, Tt, ATinfo, Ainfo, sigma, epsilon, q, args.ss, args.qs, new_top_name)
    E_SEA=ESEA(args.trj_l, top_file, args.tpr_l,new_gro_name)
    E_LIQ=EGromacs(args.trj_l, top_file, gompp_name, args.tpr_l, 'non-Water') 
    E_VAC=EGromacs(args.trj_v, top_file, gompp_name, args.tpr_v, 'System')
    np.savetxt('E_SEA',E_SEA,fmt='%10.5f')
    np.savetxt('E_LIQ',E_LIQ,fmt='%10.5f')
    np.savetxt('E_VAC',E_VAC,fmt='%10.5f')

if __name__ == '__main__': #Weird Python way to execute main()
    main()



