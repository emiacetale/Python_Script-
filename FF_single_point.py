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

top_num ='{}_{:02d}_{:02d}.top'

def read_top(topology):
    with open(topology,'r') as f:
         TOP = f.readlines()
    #print TOP
    Top=[]
    for l in TOP:
        if l[0]!=";" and l!="\n" and l[0]!="#":
            Top.append(l)
    #print Top
    for i, line in enumerate(Top):
        if 'atomtypes' in line: 
            Atb=i+1
        if 'moleculetype' in line: 
            Ate=i
        if 'atoms'     in line: 
            Ab=i+1
        if 'bonds'     in line: 
            Ae=i
            break
    #untouched pieces (head, center, tail) ot the topology 
    #print Atb,Ate,Ab,Ae
    Th=Top[:Atb]
    Tc=Top[Ate:Ab]
    Tt=Top[Ae:-1] #We remove the presence of water from the topology 

    #Things that needs to be changed (where there is sigma and q)
    At=Top[Atb:Ate]
    A=Top[Ab:Ae]
    #We can return directly the segments, but this is more readable 
    return Th, Tc, Tt, At, A

def get_s_e(At):
    info=[]
    sigma=[] 
    epsilon=[] 
    for line in At:
        line=line.split()
        #print line
        I='{:3s}  {:3s}   {:9.6f}  {:9.6f}  {:1s} '.format(line[0],line[1],float(line[2]),float(line[3]),line[4])
        info.append(I)
        sigma.append(float(line[5]))
        epsilon.append(float(line[6]))
    #print sigma
    #print epsilon
    return info, sigma, epsilon

def get_q(A):
    info=[]
    q=[]
    for line in A:
        line=line.split()
        q.append(float(line[6]))
        I='  {:3d}    {:>3s}  {:3d}   {:>3s}   {:>3s}  {:2d}   '.format(int(line[0]),line[1],int(line[2]),line[3],line[4],int(line[5]))
        info.append([I,float(line[7])])
    return info, q

def write_top(Th, Tc, Tt, ATinfo, Ainfo, s, e, q, ss, qs, filename):
    with open(filename,'w') as outfile:
         for line in Th: 
             print >>outfile, line[:-1] 
         for at_i, sigma, epsilon in zip(ATinfo, s, e):
             line=at_i+'{:12.9f}    {:12.9f} '.format((sigma+ss*sigma),epsilon)
             print >>outfile, line
         for line in Tc:
             print >>outfile, line[:-1]
         for a_i, Q in zip(Ainfo, q):
             line='{:8.5f}   {:8.5f}'.format((Q+qs*Q),a_i[1])
             line=a_i[0]+line
             print >>outfile, line
             #print >>outfile, "\n"
         for line in Tt:
             print >>outfile, line[:-1]
         print >>outfile, "\n"
    return [filename]

def create_gro(name, group_out, trj_name, tpr_name):
    p=sp.Popen('trjconv -f '+trj_name+' -o '+name+'.gro -s '+tpr_name+' -sep -pbc mol 2>/dev/null', shell=True, stdin=sp.PIPE)
    p.stdin.write(group_out+' \n')
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('ls '+name+'*.gro',shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    return output.split() 

def energy_SEA(filename_g):
    sp.call(["cp",filename_g,"gromacs.gro"])
    #p=sp.Popen("/home/ebrini/software/SEA_LIBO/bin/solvate -s gromacs -ce none -d 12 -i 500 2>/dev/null", shell=True, stdout=sp.PIPE)
    #p.communicate()
    p=sp.Popen("/home/ebrini/software/SEA/bin/solvate -s gromacs -ce none -d 12 -i 500", shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    output=output.split('\n')
    for line in output:
        #print line
        #if 'Non-Polar' in line: 
        if 'Total' in line: 
            #print line
            line=line.split()
            Enp=float(line[1])
    #print Enp
    #p=sp.Popen("/home/ebrini/software/FSEA/FSEA_adp_big_mol.exe < surface.povdat", shell=True, stdout=sp.PIPE)
    #(output, err) =  p.communicate()
    #print Enp
    #print output 
    #print output
    #Ep=float(output)
    #p=sp.Popen("rm surface.povdat", shell=True)
    #p.communicate()
    return Enp #((Ep/4.184)+Enp) #Our UOM is kcal wile Libo's FSEA is in KJ

def Eone(gro_files, top_file):
    sp.call(["cp",top_file,"gromacs.top"])
    E=[]
    for gro_f in gro_files:
         E.append(energy_SEA(gro_f))
    return np.array(E)

def Eall(gro_files, top_files):
    E=[]
    for top_file in top_files:
        print top_file
        e=Eone(gro_files, top_file)
    	E.append(e)
    return E 

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-trj_l', type=str, help='gromacs trajectory of the solute in solution to read')
    parser.add_argument('-trj_v', type=str, help='gromacs trajectory of the solute in vacuum to read')
    parser.add_argument('-top', type=str, help='gromacs topology to read')
    parser.add_argument('-tpr', type=str, help='gromacs tpr to read')
    parser.add_argument('-qs', type=float, help='scaling of the charge with respect to the original FF')
    parser.add_argument('-ss', type=float, help='scaling of the sigma  with respect to the original FF')
    return parser.parse_args()

def main():
    #"Hard-coded" variables: not much sense to change them at this point, but it is nice to have an handle for future 
    new_top_name="TOP.top"     # Name of topology file 
    new_gro_name="struct"  # Name of conf file 
    kBT=300*1.9872041E-3   # KT in kcal/mol

    #Script:
    args = parse_args()                                                                      #We read some input 
    Th, Tc, Tt, At, A=read_top(args.top)                                                        #Read and decompose the 
    ATinfo, sigma, epsilon  = get_s_e(At)                                                    #    topology
    Ainfo, q=get_q(A)                                                                        #    ...
    top_file=write_top(Th, Tc, Tt, ATinfo, Ainfo, sigma, epsilon, q, args.ss, args.qs, new_top_name)
    gro_files=create_gro(new_gro_name, 'non-Water', args.trj_l, args.tpr) #Create gro files

    #E0=Eone(gro_files, top_num.format(new_top_name, zero[0], zero[1])) # E of each config for original FF 
    #Ea=Eall(gro_files, top_files)                        # " "  "    "      "   all  
    #dg0=DG0(E0,nbin)                                     # DG* original FF
    #ddg_all= DDG(Ea,E0,kBT)                              # DG* all 
    #dg_all=ddg_all+dg0                                   # calc DG*-DG*_target 

if __name__ == '__main__': #Weird Python way to execute main()
    main()


