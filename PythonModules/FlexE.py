import numpy as np
import mdtraj as md 
from tqdm import tqdm as tqdm
from joblib import Parallel, delayed
from copy import deepcopy as dc

#Function to compute FlexE over all the frames in a trajectory given a reference structure
#Trajectory and reference should have the same number of atoms and only the important one (eg CA)
#The paramenters are:
# trajectory -positional argument-
# reference structure -positional argument-
# k_b = spring constant for bonded interactions in [kcal / (mol A^2)] 
# k_nb = spring k for non-bonded int in [kcal / (mol A^2)] after divinding for a distance and elevating to the 6th
# scaling_factor = scaling factor for spring constants (0.4 in FlexE paper)
# ij_bond = max bond separation -contact order- to use k_b instead of k_nb
# r_nb_max = only residues closer than r_nb_max [A] will be considered 
# pairs_to_exclude = interaction to exclude as pair of atoms [[i,j],[k,l],...] !! NB: i < j and k < l !!
#

#Computes the FlexE for all frames of a trajectory vs a reference
def FlexE_1D(traj,ref,k_b=60,k_nb=6,scaling_factor=0.4,ij_bond=3,r_nb_max=12,pairs_to_exclude=[]):
    traj=dc(traj)
    ref=dc(ref)
    EB=E_bond(traj,ref,k_b,scaling_factor,ij_bond,pairs_to_exclude)
    ENB=E_nonbond(traj,ref,k_nb,scaling_factor,ij_bond,r_nb_max,pairs_to_exclude)
    E=(EB+ENB)/ref.top.n_residues
    del traj
    del ref
    return(np.array(E))

#This function is here for convenience to compute the 2D FlexE (all frames vs each others) in parallel
def FlexE_2D(traj,k_b=60,k_nb=6,scaling_factor=0.4,ij_bond=3,r_nb_max=12.0,pairs_to_exclude=[],nj=8):
    frame_inexes=np.arange(0,len(traj))
    RMSD2D=Parallel(n_jobs=nj)(delayed(FlexE_1D)(traj,
                                                 traj[i],
                                                 k_b=k_b,
                                                 k_nb=k_nb,
                                                 scaling_factor=scaling_factor,
                                                 ij_bond=ij_bond,
                                                 r_nb_max=r_nb_max,
                                                 pairs_to_exclude=pairs_to_exclude) for i in tqdm(frame_inexes))
    return(np.array(RMSD2D))

######################################################
# Below are functions called by FlexE_1D             #
######################################################


#FlexE contribution for residues closer than 3 bonds
def E_bond(traj,r,k_b,scaling_factor,ij_bond,pairs_to_exclude):
    pairs=get_pair_bonded(r.top.n_residues,ij_bond)
    if len(pairs_to_exclude)>0:
        pairs=exclude_pairs(pairs,pairs_to_exclude)
    springs=[scaling_factor*(k_b/((p[1]-p[0])*(p[1]-p[0]))) for p in pairs]
    E=get_E(traj,r,pairs,springs)
    del pairs
    del springs
    return(E)

#FlexE contribution for residues further than 3 bonds
def E_nonbond(traj,r,k_nb,scaling_factor,ij_bond,r_nb_max,pairs_to_exclude):
    pairs=get_pair_nb(r,ij_bond,r_nb_max)
    if len(pairs_to_exclude)>0:
        pairs=exclude_pairs(pairs,pairs_to_exclude)
    springs=get_springs_nb(r,pairs,k_nb)*scaling_factor
    E=get_E(traj,r,pairs,springs)
    del pairs
    del springs
    return(E)

#This energy function computes k_springs*(dist^2). It is the same for bonded and non bonded 
def get_E(traj,r,pairs,springs):
    d0=get_dist_A(r,pairs)[0]
    dt=get_dist_A(traj,pairs)
    dt=np.square(dt-d0)
    dt=dt*springs
    E=np.sum(dt,axis=1)
    del d0
    del dt
    return(0.5*E) #The 0.5 accounts for the 1/2 factor in the spring equation

#Bonded residues are defined as residues with a contact order (CO) <= ij_bond 
def get_pair_bonded(nr,ij_bond):
    pb=[]
    for i in np.arange(0,nr):
        for j in np.arange(i+1,i+ij_bond+1):
            if j<nr:
                pb.append([i,j])
    return(pb)

#NB-residue pairs are all residues that have CO>ij_bond and that are inside a cut off r_nb_max
def get_pair_nb(r,ij_bond,r_nb_max):
    all_pairs=r.top.select_pairs(np.arange(0,r.top.n_atoms),np.arange(0,r.top.n_atoms))
    dist=get_dist_A(r,all_pairs)[0]
    pairs=[ p for p,d in zip(all_pairs,dist) if (d<r_nb_max and (p[1]-p[0])>ij_bond)]
    del all_pairs
    del dist
    return(pairs)

#Sometimes we want to exclude residue pairs from the FlexE calculation (e.g. bonded interactions between monomers, ...
# ... or bonds in drugs). This is a flexible approach (but not super fast if the list is long)
def exclude_pairs(pairs,pairs_to_exclude):
    P=np.array(pairs)
    pte=np.array(pairs_to_exclude)
    ok_pairs=[]
    for p in P:
        if not(np.any([np.all(t==p) for t in pte])):
            ok_pairs.append(p)
        #else:
        #    print("excluding: ",p)
    del P
    del pte
    return(np.array(ok_pairs))
    
#NB-springs are function of the residue distance in the referece
def get_springs_nb(r,pairs,k_nb):
    tmp=(np.ones(len(pairs))*k_nb)/get_dist_A(r,pairs)[0]
    return(np.power(tmp,6.0*np.ones(len(pairs))))
    #return(np.square((np.ones(len(pairs))*k_nb)/get_dist_A(r,pairs)[0]))

#mdtraj gives distances in nm, we have this function just to get them in A 
def get_dist_A(t,p,opt=True,per=False):
    return(md.compute_distances(traj=t,atom_pairs=p,opt=opt,periodic=per)*10.0)


