#This module defines few useful functions to run SEA calulations in a python script:
#
# create_gro --> wrap a trjconv command to create the grefile of the solute from a gromacs traj
# clean_topology --> gives back a topology file without references to water
# energy_SEA --> gives back sea solvation free energy given a gro and a top file 
#
#

#import mdtraj as md
import subprocess as sp
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

def energy_SEA(filename_g,filename_top,quadrupole=False,fieldsea=False,detail_np=12, 
           iteration_sea= 50, iteration_fsea=1, surfcull_fsea=5, detail_fsea=24, 
           remove_gro_file=True):
    solvate="/home/ebrini/software/newSEA/bin/solvate "
    if remove_gro_file:
        sp.call(["mv",filename_g,"gromacs.gro"])
    else:
        sp.call(["cp",filename_g,"gromacs.gro"])
    #We cp the top file to a standard name 
    sp.call(["cp",filename_top,"gromacs.top"])

    #Dipolar SEA
    if (not(quadrupole) and not(fieldsea)):
        line=solvate+"-s gromacs -ce rfgb -d "+str(detail_np)+" -i "+str(iteration_sea)+" "
    #Quadrupole SEA
    if quadrupole:
      line=solvate+"-s gromacs -ce rfgb -q -d "+str(detail_np)+" -i "+str(iteration_sea)+" "
    #Field SEA (NP only)
    if fieldsea:
      line=solvate+"-s gromacs -f -d "+str(detail_np)+" -i "+str(iteration_sea)+ \
           " -fd "+str(detail_fsea)+" -fbi "+str(surfcull_fsea)+" -fi "+str(iteration_fsea)+" "  
    p=sp.Popen(line, shell=True, stdout=sp.PIPE)
    (output, err) =  p.communicate()
    output=output.split('\n')
    for line in output:
        if 'Non-Polar' in line:
            line=line.split()
            Enp=float(line[1])
        if 'Total' in line: 
            line=line.split()
            Etot=float(line[1])
    return Etot,Enp,Etot-Enp #We return DG*, DG*_np, DG*_p

def clean_topology(top,name_new_top="top_clean.top"):
    new_top=[]
    with open(top,'r') as F_in: 
      with open(name_new_top,'w') as F_out:
        for line in F_in: 
          if not(('ions.itp' in line) or ('posre.itp' in line) or ('tip3p.itp' in line) or ('SOL' in line)):
             F_out.write(line)
    return name_new_top

if __name__ == '__main__': #Weird Python way to execute main()
    main()



