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
import cPickle as pic
import os 

def create_gro(name, group_out, trj_name, tpr_name):
    p=sp.Popen('trjconv -f '+trj_name+' -o '+name+'.gro -s '+tpr_name+' -sep -pbc mol &>/dev/null', shell=True, stdin=sp.PIPE)
    p.stdin.write(group_out+' \n')
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('ls '+name+'*.gro',shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    list=output.split() 
    return [name+str(i)+'.gro' for i in np.arange(0,len(list),1,dtype=int)] #Orrible but effective way to create an ordered

def energy_SEA(filename_g,filename_top,quadrupole):
    sp.call(["cp",filename_g,"gromacs.gro"])
    sp.call(["cp",filename_top,"gromacs.top"])
    if quadrupole:   
      line="/home/ebrini/software/SEA/bin/solvate -s gromacs -ce none -d 12 -i 500  -q #2> /dev/null"
    else:
      line="/home/ebrini/software/SEA/bin/solvate -s gromacs -ce none -d 12 -i 500 -q #2> /dev/null"
    p=sp.Popen(line, shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    output=output.split('\n')
    for line in output:
        if 'Total' in line: 
            line=line.split()
            Etot=float(line[1])
        if 'Non-Polar' in line:
            line=line.split()
            Enp=float(line[1])
    return Etot,Enp,Etot-Enp 

def energy_fSEA(filename_g, filename_top):                    
    sp.call(["cp",filename_g,"gromacs.gro"])
    sp.call(["cp",filename_top,"gromacs.top"])
    p=sp.Popen("/home/ebrini/software/SEA_LIBO/bin/solvate -s gromacs -ce none 2>/dev/null", shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    output=output.split('\n')
    for line in output:
        if 'Non-Polar' in line: 
            line=line.split()
            Enp=float(line[1])

    p=sp.Popen("/home/ebrini/software/FSEA/FSEA_adp_big_mol.exe < surface.povdat", shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()           
    #print output                              
    Ep=float(output)                           
    p=sp.Popen("rm surface.povdat", shell=True)
    p.communicate()                            
    Ep=Ep/4.184 #Our UOM is kcal wile Libo's FSEA is in KJ
    return Ep+Enp,Enp,Ep

def clean_topology(top):
    name_new_top='top_clean.top'
    new_top=[]
    with open(top,'r') as F_in: 
      with open(name_new_top,'w') as F_out:
        for line in F_in: 
          if not(('ions.itp' in line) or ('posre.itp' in line) or ('tip3p.itp' in line) or ('SOL' in line)):
             F_out.write(line)
    return name_new_top

def clean_folder(file_list):
    for file in file_list: 
        os.remove(file)


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-trj', type=str, help='gromacs trajectory of the solute in solution to read')
    parser.add_argument('-top', type=str, help='gromacs topology to read (without water!)')
    parser.add_argument('-tpr', type=str, help='gromacs tpr to read for solution sim')
    parser.add_argument('-fs', help='perform a field-sea calc instead of dipolar-SEA',
                              action='store_true', default=False)
    parser.add_argument('-qs', help='perform a quadrupole-sea calc instead of dipolar-SEA',
                              action='store_true', default=False)

    return parser.parse_args()

def main():
    new_gro_name='struct'
    args = parse_args() #We read some input 

    #Script:
    if args.fs and args.qs:
        print '!!! Error !!!'
        print 'You may request only one type of SEA calc at a time'
        exit()
   
    if args.fs:   #We need to copy some files to have FSEA working
       p=sp.Popen("cp /home/ebrini/software/FSEA/water_dist.csv /home/ebrini/software/FSEA/pos_formula.dat /home/ebrini/software/FSEA/neg_formula.dat .", shell=True, stdout=sp.PIPE)

    E=[]
    Enp=[]
    Ep=[]
    #print args.trj[-3:]
    #exit() 

    if args.trj[-3:]=='trr' or args.trj[-3:]=='trj':
        gro_files=create_gro(new_gro_name, 'non-Water', args.trj, args.tpr) #Extract traj in gro files 
    else: 
        gro_files=[args.trj]
    
    cleaned_top=clean_topology(args.top)
    
    for coord_file in  gro_files:
      #print coord_file
      if args.fs or args.qs:  #we want to calulate quadrupole or field sea
          if args.fs:         #Field SEA
                tmp=energy_fSEA(coord_file,cleaned_top)
                
                #print "FSEA"
          else:               #Quadrupole SEA
                tmp=energy_SEA(coord_file,cleaned_top,True)
                #print "QSEA"
      else:                   #Dipole SEA
                tmp=energy_SEA(coord_file,cleaned_top,False)
                #print "DSEA"
      E.append(tmp[0])
      Enp.append(tmp[1])
      Ep.append(tmp[2])

    E=np.array(E)
    Enp=np.array(Enp)
    Ep=np.array(Ep)
    
    if args.fs or args.qs:  #Make a BU file with the enerdy and an appropriate name
       if args.fs:         #Field SEA
          with open('FSEA_E.pcl','wb') as fd: pic.dump((E,Enp,Ep),fd) 
       else:               #Quadrupole SEA
          with open('QuadSEA.pcl','wb') as fd: pic.dump((E,Enp,Ep),fd) 
    else:                   #Dipole SEA
          with open('DipSEA.pcl','wb') as fd: pic.dump((E,Enp,Ep),fd) 

    print 'DG*   = '+str(np.average(E))+' kcal/mol'
    print 'DGnp* = '+str(np.average(Enp))+' kcal/mol'
    print 'DGp*  = '+str(np.average(Ep))+' kcal/mol'

    if args.trj[-3:]!='gro':
        clean_folder(gro_files)

if __name__ == '__main__': #Weird Python way to execute main()
    main()



