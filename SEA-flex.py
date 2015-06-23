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

def create_mdp(name):
    A=[';MD param','integrator=md','dt=0.002 ;','nsteps=2500000 ; 5 ns',';COM motion removal','comm-mode=linear','nstcomm=100',
       ';Output','nstxout=0','nstvout=0','nstfout=0','nstcalcenergy=1','nstenergy=1','nstlog=0',';pbc','pbc=xyz',
       ';Electrostatic','coulombtype=cutoff','rcoulomb=1.5','vdwtype=cut-off','rvdw=1.5','rlist=1.5',';dispcorr=EnerPres',
       ';T copupling','tcoupl=no','pcoupl=no','compressibility=4.5E-5','constraints=h-bonds']
    with open(name,'w') as f:
         for line in A:
             print >>f, line

def create_gro(name, group_out, trj_name, tpr_name):
    p=sp.Popen('trjconv -f '+trj_name+' -o '+name+'.gro -s '+tpr_name+' -sep -pbc mol &>>log', shell=True, stdin=sp.PIPE)
    p.stdin.write(group_out+' \n')
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('ls '+name+'*.gro',shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    list=output.split() 
    return [name+str(i)+'.gro' for i in np.arange(0,len(list),1,dtype=int)] #Orrible but effective way to create an ordered

def ESEA(traj, top, tpr, new_gro_name):
    gro_files=create_gro(new_gro_name, 'Protein', traj, tpr) #Create gro files
    E=Eone(gro_files, top)
    return (E) 

def energy_SEA(filename_g):
    sp.call(["mv",filename_g,"gromacs.gro"])
    p=sp.Popen("/home/ebrini/software/SEA/bin/solvate -s gromacs -ce none -d 12 -i 500 2>>log", shell=True, stdout=sp.PIPE)
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
    p=sp.Popen('trjconv -f '+traj+' -o traj0.trr -s '+tpr+' -pbc mol &>>log ', shell=True, stdin=sp.PIPE)
    p.stdin.write(group+' \n') #We should be sure to remove water from the traj
    p.communicate()[0]            
    p.stdin.close()
    p=sp.Popen('grompp -f '+grompp+' -p '+top+' -c gromacs.gro &>>log ', shell=True, stdin=sp.PIPE)
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('mdrun -rerun traj0.trr -nt 1 &>>log ', shell=True, stdin=sp.PIPE)
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('g_energy -xvg none -o energy.xvg &>>log ', shell=True, stdin=sp.PIPE)
    p.stdin.write('Potential \n')
    p.stdin.write('\n')
    p.communicate()[0]
    p.stdin.close()
    X, E=np.loadtxt('energy.xvg', dtype=float, unpack=True)
    p=sp.Popen('rm energy.xvg ener.edr md.log mdout.mdp topol.tpr traj.trr traj0.trr', shell=True, stdin=sp.PIPE)
    p.communicate()[0]
    p.stdin.close()
    return (E*0.239005736)                     #We convert kJ in kcal

def calc_DE(SEA_l,SEA_v,E_l,E_v): 
    L=SEA_l+E_l
    V=SEA_v+E_v
    de=[ np.average(L-ev) for ev in V ]
    return de

def plot_distr(E, kT, nsteps):
    prob,edges=np.histogram(E,bins=nsteps,density=True)
    step=edges[1]-edges[0]
    edges=edges[:-1]-0.5*step
    
    fig=plt.figure()
    ax1 = fig.add_axes([0.25, 0.15, 0.65, 0.8])
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel(r'$\Delta G^* / kcal \, mol^{-1}$',fontsize=20)
    ax1.set_ylabel(r'$P(\Delta G*$)',fontsize=20)
    p1,=ax1.plot(edges,prob, color='blue', linewidth=2, linestyle="-", label=r'$P(\Delta G* )')

    prob=[ p*np.exp(-e/kT) for p,e in zip(prob,edges) ]
    prob=prob/np.trapz(prob,edges)

    p2,=ax1.plot(edges,prob, color='orange', linewidth=2, linestyle="-", label=r'$e^{- \beta \Delta G} P(\Delta G* )')

    fig.savefig('P_distr.png') 




def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-trj_l', type=str, help='gromacs trajectory of the solute in solution to read')
    parser.add_argument('-trj_v', type=str, help='gromacs trajectory of the solute in vacuum to read')
    parser.add_argument('-top', type=str, help='gromacs topology')
    parser.add_argument('-tpr_l', type=str, help='gromacs tpr to read for solution sim')
    parser.add_argument('-tpr_v', type=str, help='gromacs tpr to read for vac sim')
    parser.add_argument('-pcl', help='load E from a pickle file instead of computing them',
                         action='store_false',default=True)
    parser.add_argument('-pcl_file', type=str, help='pickle-file to load/save',default='E_data_BU') 
    return parser.parse_args()

def main():
    #"Hard-coded" variables: not much sense to change them at this point, but it is nice to have an handle for future 
    new_gro_name="struct"  # Name of conf file 
    gompp_name="grompp.mdp"
    kT=0.59616123 #kT in kcal/mol


    #Script:
    args = parse_args()                                                                      #We read some input
    create_mdp(gompp_name) 
    if args.pcl:
       SEA_l=ESEA(args.trj_l, args.top, args.tpr_l,new_gro_name)
       print "SEA_l"
       SEA_v=ESEA(args.trj_v, args.top, args.tpr_v,new_gro_name)
       print "SEA_v"
       E_l=EGromacs(args.trj_l, args.top, gompp_name, args.tpr_l, 'non-Water') 
       print "E_l"
       E_v=EGromacs(args.trj_v, args.top, gompp_name, args.tpr_v, 'System')
       print "E_v"
       with open(args.pcl_file,'wb') as fd: pic.dump((SEA_l,SEA_v,E_l,E_v),fd)
    else:
       with open(args.pcl_file,'rb') as fd: SEA_l,SEA_v,E_l,E_v=pic.load(fd)
    
    temp=[ np.exp(-E/kT) for E in SEA_v ]
    DG=-kT*np.log(np.average(temp))
    print "DG     = "+str(DG)+" kcal/mol"
    plot_distr(SEA_v, kT, 100)


if __name__ == '__main__': #Weird Python way to execute main()
    main()



