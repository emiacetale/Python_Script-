#!/usr/bin/env python

#import mdtraj as md
import numpy as np
import subprocess as sp

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

def get_s_e(At):
    info=[]
    sigma=[] 
    epsilon=[] 
    for line in At:
        line=line.split()
	I='{:3s}  {:3d}   {:9.6f}  {:9.6f}  {:1s} '.format(line[0],int(line[1]),float(line[2]),float(line[3]),line[4])
	info.append(I)
	sigma.append(float(line[5]))
	epsilon.append(float(line[6]))
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

def create_gro(name, group_out):
    p=sp.Popen('trjconv -f traj.trr -o '+name+'.gro -s topol.tpr -sep -e 200', shell=True, stdin=sp.PIPE)
    p.stdin.write(group_out+' \n')
    p.communicate()[0]
    p.stdin.close()
    p=sp.Popen('ls '+name+'*.gro',shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    return output.split() 

def Eone(gro_files, top_file):
    sp.call(["cp",top_file,"gromacs.top"])
    E=[]
    for gro_f in gro_files:
        sp.call(["cp",gro_f,"gromacs.gro"])
	p=sp.Popen("/home/ebrini/software/SEA/bin/solvate -s gromacs -d 12 -i 500 2>/dev/null", shell=True, stdout=sp.PIPE)
	(output, err) =  p.communicate()
        for line in output.split('\n'):
	  if "Total" in line:
	     line=line.split()
	     Etot=float(line[1])
	     E.append(Etot)
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


def main():
    #Variables: 
    qrange=[-0.1,0.1]
    srange=[-0.1,0.1]
    #npoints=21 #MUST be odd!!!!
    npoints=3
    new_top_name="TOP"
    new_gro_name="struct"
    nbin=50
    kBT=300*1.9872041E-3 #in kcal/mol
    dg_targhet=-4.756 #kcal

  
    #Script:
    qrange=np.linspace(qrange[0],qrange[1],num=npoints)
    srange=np.linspace(srange[0],srange[1],num=npoints)
    Th, Tc, Tt, At, A=read_top("topol.top")
    ATinfo, sigma, epsilon  = get_s_e(At)
    Ainfo, q=get_q(A)
    gro_files=create_gro(new_gro_name, 'non-Water')
    top_files=write_tops(Th, Tc, Tt, ATinfo, Ainfo, sigma, epsilon, q, srange, qrange, new_top_name, npoints)
    E0=Eone(gro_files, top_files[int(len(top_files)/2)]) #We pass all the config file and only the middle top (original FF)
    Ea=Eall(gro_files, top_files) 
    dg0=DG0(E0,nbin)
    ddg_all= DDG(Ea,E0,kBT)
    dg_all=ddg_all+dg0
    pint_res(dg_all, qrange, srange, dg_targhet)
    #print dg_all


    #print dg0

if __name__ == '__main__': #Weird Python way to execute main()
    main()



