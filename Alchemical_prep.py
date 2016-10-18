#!/usr/bin/env python

#import mdtraj as md
import argparse
import numpy as np
import subprocess as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

top_num ='{}_{:02d}_{:02d}.top'


def read_top(topology):
    with open(topology,'r') as f:
         TOP = f.readlines()
    #print TOP
    Top=[]
    for l in TOP:
        if l[0]!=";" and l!="\n":
    	   Top.append(l)
    for i, line in enumerate(Top):
        if 'atomtypes' in line: 
	    Atb=i+1
	if 'pairtypes' in line: 
	    Ate=i
        if 'atoms'     in line: 
	    Ab=i+1
	if 'bonds'     in line: 
	    Ae=i
	    break
    #untouched pieces (head, center, tail) ot the topology 
    Th=Top[:Atb]
    Tc=Top[Ate:Ab]
    Tt=Top[Ae:-1] #We remove the presence of water from the topology 

    #Things that needs to be changed (where there is sigma and q)
    At=Top[Atb:Ate]
    A=Top[Ab:Ae]
    #We can return directly the segments, but this is more readable 
    return Th, Tc, Tt, At, A


#	print >>out ,'{:4.2f} {:4.2f} {:9.6f} {:9.6f}'.format( q, s, dg, dg-dgt)   

def read_pdb(pdbf,nres):
    with open(pdbf) as F:
         temp=F.readlines()
    b=100000
    e=0
    for i, line in enumerate(temp):
        if line[0:4]=='ATOM':
           line=line.split()
           if line[5]==nres:
              if i<b: b=i
              if i>e: e=i
    e=e+1
    res=[line.split() for line in temp[b:e]]
    return temp[:b],res,temp[:e]

def read_rtp(rtpf,rn1,rn2):
    with open(rtpf,'r') as F:
         temp=F.readlines()
    breaks=[]                         #We have to find the beginning and the end of 
    for i,line in enumerate(temp):    #each AA. Therefore I am looking for all the lines
        if line[0]=='[':              #that start with a square bracket...
            breaks.append(i)
    b1=0
    e1=0
    for i,b in enumerate(breaks):     #...Then I look for a specific residue name in each  
        if '[ '+rn1 in temp[b]:       #of these lines and I extract the section begin
            b1=breaks[i]              #and end line...
            e1=breaks[i+1]
    b2=0
    e2=0
    for i,b in enumerate(breaks):
        if '[ '+rn2 in temp[b]:
            b2=breaks[i]
            e2=breaks[i+1]

    rtp1=temp[b1:e1]                  #...And I use these to pass the tpr of the 2 residue               
    rtp2=temp[b2:e2]
    return temp, rtp1, rtp2

def get_rtp_dict(r1):
    breaks=[]
    for i,line in enumerate(r1):
        if '[' in line:
            breaks.append(i)
    d={}
    for i,b in enumerate(breaks[:-1]):
        d[r1[b][3:-3]]=r1[b+1:breaks[i+1]]
    d[r1[breaks[-1]][3:-3]]=r1[breaks[-1]+1:]
    return d

def rename_atom(s):
    if s[-1].isalpha():
        return s+'K'
    else:
        return s[:-1]+str(int(s[-1])+5)

def rename_rtp_atom(lines):
    #print lines
    atoms=[]
    for line in lines:
        line=line.split()
        line[0]=rename_atom(line[0])
        atoms.append(["   "+line[0]+"   "+line[1]+"   "+line[2]+"   "+line[3]+"   "])
    print atoms
    return atoms


def combine_rtp(r1,r2):
    d_r1=get_rtp_dict(r1)
    d_r2=get_rtp_dict(r2) 
    #for key,value in d_r1.iteritems():
    #    print key
    #for key,value in d_r2.iteritems():
    #    print key
    d_r2['atoms']=rename_rtp_atom(d_r2['atoms'])
    return 1 


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb1', help='1st PDB file')
    parser.add_argument('-pdb2', help='2nd PDB file')
    parser.add_argument('-rtp', help='residue topology file')
    parser.add_argument('-nmut', help='residue to mutade (PDB numbering)')
    return parser.parse_args()

def main():
    args = parse_args() #We read some input 
    head1,res1,tail1=read_pdb(args.pdb1,args.nmut)
    head2,res2,tail2=read_pdb(args.pdb2,args.nmut)
    resname1=res1[0][3]
    resname2=res2[0][3]
    rtp_old,rtp1,rtp2=read_rtp(args.rtp, resname1, resname2)
    rtp_new=combine_rtp(rtp1[1:],rtp2[1:])       #We remove the name of the residue 




if __name__ == '__main__': #Weird Python way to execute main()
    main()



