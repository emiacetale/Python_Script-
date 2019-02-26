#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import cPickle as pickl
from joblib import Parallel, delayed
import copy as copy 

def calc_RMDS(trj,ref):
    coord_trj=trj.xyz
    coord_ref=ref.xyz[0]
    dist=np.square(coord_trj-coord_ref)
    RMSDs=np.array([np.sqrt((np.sum(d))/len(d)) for d in dist])
    return RMSDs
       
def wrap_RMSD(T, frames, P_index, L_index):
    RMSDs=[]
    t=copy.deepcopy(T) # we need to deep-copy the traj to edit it in parallel
    for i in frames:
        print i 
        r=t[i]
        t=t.superpose(r, frame=0, atom_indices=P_index, 
                   ref_atom_indices=P_index, parallel=False)
        RMSDs.append(calc_RMDS(t.atom_slice(L_index),r.atom_slice(L_index)))
    return np.array(RMSDs)

def load_trj(trj_path,stride,start,top):
    if top==None:                        #If no top is given 
        t=md.load(trj_path)
    else:                                #if top is given
        t=md.load(trj_path,top=top)
    #print(t)
    if t.n_frames < start:
        print '!!!ERROR!!! There are NOT enough frames in the traj' 
        print 'There are only {:} frames in this trajectory'.format(t.n_frames)
        print '!!ABORTING!!!'
        exit()
    return t.atom_slice(t.top.select('name CA'))[start::stride]

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-n_proc',type=int, help='number of parallel processes you spawn', required=True)
    parser.add_argument('-i_proc',type=int, help='which process is this one (1 based index)', required=True)
    parser.add_argument('-n_cores',type=int, help='number of cores available to each processes', required=True)
    parser.add_argument('-trjs_name',type=str, help='name of the trajectory [trajectory]', default='trajectory')
    parser.add_argument('-top_file',type=str, help='name of the topology file (if needed, eg for dcd traj) [None]')
    parser.add_argument('-trjs_ext',type=str, help='trajectory file-type (extension) [h5]', default='h5')
    parser.add_argument('-n_trajs',type=int, help='numer of trajectoryes to load [5]', default=5)
    parser.add_argument('-trj_eq_steps',type=int, help='initial number of steps to ignore from the trajectory [5000]', default=5000)
    parser.add_argument('-trj_stride',type=int, help='read every trj_stride steps [10]', default=10)
    parser.add_argument('-p1_ai',type=int, help='If passed this is the index of the CA used for alignong the structure (0 based) [None == P1]', nargs = "+")
    parser.add_argument('-p2_ai',type=int, help='If passed this is the index of the CA used to compute the RMSD (0 based) [None == P2]', nargs = "+")

    return parser.parse_args()

def main():
    args = parse_args()
    if args.n_trajs>args.n_cores:
        print '!!!ERROR!!!  Not enough cores to read all the trajs'
        print '{:} cores < {:} trajs'.format(args.n_cores,args.n_trajs)
        print '!!ABORTING!!!'
        exit()
    trajs=Parallel(n_jobs=args.n_trajs)(delayed(load_trj)('{:}.{:02d}.{:}'.format(args.trjs_name,i,args.trjs_ext),args.trj_stride,args.trj_eq_steps,args.top_file) for i in range(0,args.n_trajs))
    trajs=md.join(trajs)

    i_proc=args.i_proc-1                                                      #From 1 based (bash) to 0 based (python) index
    sep=np.linspace(0,trajs.n_frames,args.n_proc+1,endpoint=True, dtype=int)  #We divide the frames between the different parallel process so ...
    begin=sep[i_proc]                                                         # ... we know where our trajectory slice begin ...
    end=sep[i_proc+1]                                                         # ... and ends. 
    frames=np.arange(begin,end)                                               #We can now have a list of all our ref frames ...
    frames=np.array_split(frames,args.n_cores)                                # ... and split it between the cores on this node
    
    if args.p1_ai == None:                                                    #We need the atom index, that we could have passed or not  
        ai_0=trajs.top.select('chainid 0')                                    #ai_0 are used to align
    else:                                                                     #ai_1 are used to compute RMSD 
        ai_0=np.array(args.p1_ai)

    if args.p2_ai == None:
        ai_1=trajs.top.select('chainid 1')
    else:
        ai_1=np.array(args.p2_ai)


    #We compute the RMSD for the entire trajectory using our reference frames in parallel  
    RMSDs=Parallel(n_jobs=args.n_cores)(delayed(wrap_RMSD)(trajs,frames[i],ai_0,ai_1) for i in np.arange(0,args.n_cores))
    RMSDs=np.concatenate(RMSDs) #We get rid of the structure coming from the parallel part of the job 
    
    with open('2DRMSD_{:02d}.pickl'.format(i_proc),'wb') as fout:
        pickl.dump(RMSDs,fout)
    
    print "Done"

    



if __name__ == "__main__": main()

