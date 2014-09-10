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

def write_tops(Th, Tc, Tt, ATinfo, Ainfo, s, e, q, srange, qrange, top_name, npoints):
    top_list=[]
    for iq, dq in enumerate(qrange):  #Loop over charges
     for js, ds in enumerate(srange): #Loop over sigmas 
       filename = top_num.format(top_name, iq, js)
       top_list.append(filename)      #We create a list of the topfile so later we are more flexible 
       with open(filename,'w') as outfile:
         for line in Th: 
             print >>outfile, line[:-1] 
         for at_i, sigma, epsilon in zip(ATinfo, s, e):
             line=at_i+'{:12.9f}    {:12.9f} '.format((sigma+ds*sigma),epsilon)
             print >>outfile, line
         for line in Tc:
             print >>outfile, line[:-1]
         for a_i, Q in zip(Ainfo, q):
             line='{:8.5f}   {:8.5f}'.format((Q+dq*Q),a_i[1])
             line=a_i[0]+line
             print >>outfile, line
             print >>outfile, "\n"
         for line in Tt:
             print >>outfile, line[:-1]
         print >>outfile, "\n"
    return top_list

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

def DG0(E, nbin):
    P, bins = np.histogram(E, nbin, normed=1)
    bins=bins[:-1]
    dE=(bins[1]-bins[0])
    DG=np.sum(((bins+0.5*dE)*P)*dE)
    return DG    

def DDG(Ea,E0,kBT):
    ddg=[]
    for E in Ea:
        P, bins = np.histogram((E-E0))
        de=bins[1]-bins[0]
	e=-kBT*np.log(np.sum(de*P*np.exp(-(bins[:-1]+0.5*de)/kBT)))
	ddg.append(e)
    return ddg

def pint_res(dg_all,qrange,srange, dgt):
    DG=dg_all.reshape(len(qrange),len(srange)) # { {DG_s0, DG_s1, ... }@q0 , {DG_s0, DG_s1, ... }@q1 , ... }
    with open("results",'w') as out:
      print >>out, "#q_scale  s_scale  DG*   DG*-target"
      for q, el in zip(qrange, DG):
          for s, dg in zip(srange, el):
	      print >>out ,'{:4.2f} {:4.2f} {:9.6f} {:9.6f}'.format( q, s, dg, dg-dgt)   

def plot_res(dg, qrange, srange):
    m=np.amin(dg)
    M=np.amax(dg)
    if (-m)>M: M=-m
    dg=dg.reshape(len(qrange),len(srange)) 
    dg=np.transpose(dg)
    dq=qrange[1]-qrange[0]
    qrange=qrange-0.5*dq
    qrange=np.append(qrange, qrange[-1]+dq)
    ds=srange[1]-srange[0]
    srange=srange-0.5*ds
    srange=np.append(srange, srange[-1]+ds)
    fig=plt.figure()
    ax1 = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    ax10= fig.add_axes([0.91, 0.15, 0.02, 0.75])

    #print srange
    #print qrange
    #print dg
    
    ax1.pcolor(qrange, srange, dg, vmin=-M, vmax=M, cmap='bwr')
    ax1.set_xlabel(r'$\Delta q$',fontsize=25)
    ax1.set_ylabel(r'$\Delta \sigma$',fontsize=25)
    ax1.set_xlim(srange[0],srange[-1])
    ax1.set_ylim(qrange[0],qrange[-1])
    ax1.set_xticks(np.arange(qrange[0]+0.5*dq,qrange[-1],dq))
    ax1.set_yticks(np.arange(srange[0]+0.5*ds,srange[-1],ds))
    ax1.tick_params(axis='y', which='major', labelsize=20)
    ax1.tick_params(axis='x', which='major', labelsize=20)

    cb=mpl.colorbar.ColorbarBase(ax10, cmap='bwr', orientation='vertical', norm=mpl.colors.Normalize(vmin=-M, vmax=M),ticks=[-M, -M/2, 0, M/2, M])
    plt.savefig('Result.png')

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-np', type=int, help='number of points (21)', default=21)
    parser.add_argument('-dgt',type=float,  help='target solvation free energy')
    parser.add_argument('-trj_n', type=str, help='gromacs trajectory to read')
    parser.add_argument('-top_n', type=str, help='gromacs topology to read')
    parser.add_argument('-tpr_n', type=str, help='gromacs tpr to read')
    parser.add_argument('-q_r', type=float, help='range of charge that will be sampled (+-q)',default=0.1)
    parser.add_argument('-s_r', type=float, help='range of sigma that will be sampled (+-s)', default=0.1)
    parser.add_argument('-nbin',type=int, help='number of bins to calulate P(E), P(DE) distr', default=50)
    parser.add_argument('-split', help='split the computation on 4 subsets', action="store_true")
    parser.add_argument('-sect', type=int, help='subset to be calculated', default=5)
    return parser.parse_args()

def prepare_args(args):                                               #This fuction handles the in line input 
    npoints=args.np                                                   #
    if npoints % 2 == 0:                                              #
        print "!!The number of point must be odd!!"                    
        exit()
    dg_target=args.dgt                                                #These are not super needed 
    top_n=args.top_n                                                  #
    trj_n=args.trj_n                                                  #
    tpr_n=args.tpr_n                                                  #
    nbin=args.nbin                                                    #
    sect=args.sect                                                    #
    qrange=np.linspace(-(args.q_r),(args.q_r),num=npoints)            #
    srange=np.linspace(-(args.s_r),(args.s_r),num=npoints)            #
    zero=[int(np.floor(npoints/2)),int(np.floor(npoints/2))]          # To know which of the new TOP files has the original FF
    if args.split:
       if (sect > 4) or (sect < 1):
          print "!!!I cant calcualte a sector that does not exists (there are only 4 sect: from 1 to 4)!!!"
          exit()
       elif sect==1:                                   # This is needed: we can split the calc in 4 sectors and this handle it 
            srange=srange[0:int(np.ceil(npoints/2))+1] # The 4 sectors are arragned as follows:
            qrange=qrange[0:int(np.ceil(npoints/2))+1] #                    sigma
            zero=[len(srange)-1, len(qrange)-1]        #    |-10              0              +10|
       elif sect==2:                                   #   -+-----------------+-----------------+
            srange=srange[int(np.floor(npoints/2)):]   # -10|                 |                 |
            qrange=qrange[0:int(np.ceil(npoints/2))+1] #    |     sect 1      |  sect2          |
            zero=[0, len(qrange)-1]                    #    |  s=[-10:0]      | s=[ 0:10]       |
       elif sect==3:                                   #    |  q=[-10:0]      | q=[-10:0]       |
            srange=srange[0:int(np.ceil(npoints/2))+1] #q  0+-----------------+-----------------+
            qrange=qrange[int(np.floor(npoints/2)):]   #    |                 |                 |
            zero=[len(srange)-1, 0]                    #    |     sect 3      |  sect 4         |
       elif sect==4:                                   #    |  s=[-10:0]      | s=[ 0:10]       |
            srange=srange[int(np.floor(npoints/2)):]   # +10|  q=[ 0:10]      | q=[ 0:10]       |
            qrange=qrange[int(np.floor(npoints/2)):]   #   -+-----------------+-----------------+
            zero=[0,0]                                 # Yep, some points are calculated in 2 or more sector
    zero=[zero[1],zero[0]] # The order is q scaling and scaling in the rest of the script
    #print "q: ",qrange 
    #print "s: ",srange
    #print "zero:",zero,srange[zero[0]],qrange[zero[1]] 
    return npoints, dg_target, qrange, srange, top_n, trj_n, tpr_n, nbin, zero


def main():
    #"Hard-coded" variables: not much sense to change them at this point, but it is nice to have an handle for future 
    new_top_name="TOP"     # Name of topology file 
    new_gro_name="struct"  # Name of conf file 
    kBT=300*1.9872041E-3   # KT in kcal/mol

    #Script:
    args = parse_args()                                                                      #We read some input 
    npoints, dg_target, qrange, srange, top_n, trj_n, tpr_n, nbin, zero = prepare_args(args) #And we prepare it  
    Th, Tc, Tt, At, A=read_top(top_n)                                                        #Read and decompose the 
    ATinfo, sigma, epsilon  = get_s_e(At)                                                    #    topology
    Ainfo, q=get_q(A)                                                                        #    ...
    gro_files=create_gro(new_gro_name, 'non-Water', trj_n, tpr_n)                            #Create gro files
    top_files=write_tops(Th, Tc, Tt, ATinfo, Ainfo, sigma, epsilon, q, srange, qrange, new_top_name, npoints) #and top
    E0=Eone(gro_files, top_num.format(new_top_name, zero[0], zero[1])) # E of each config for original FF 
    Ea=Eall(gro_files, top_files)                        # " "  "    "      "   all  
    dg0=DG0(E0,nbin)                                     # DG* original FF
    ddg_all= DDG(Ea,E0,kBT)                              # DG* all 
    dg_all=ddg_all+dg0                                   # calc DG*-DG*_target 
    pint_res(dg_all, qrange, srange, dg_target)          # Print results
    plot_res(dg_all-dg_target, qrange, srange)           # Plot results

if __name__ == '__main__': #Weird Python way to execute main()
    main()



