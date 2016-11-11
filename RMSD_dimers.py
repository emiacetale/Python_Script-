#!/usr/bin/env python  

import numpy as np
import copy as copy
import argparse
import mdtraj as md 
import cPickle as pickle 

# f.write('{:5s}   {:8.6f}\n'.format(el[0],el[1]))

def get_trjs(trj_files,top='XXX'):
    trajs = []
    for el in trj_files:
        if top != 'XXX':
           tmp=md.load(el,top=top)
        else:
           tmp=md.load(el)
        #We keep only BB coordinates
        ai=tmp.top.select("backbone")
        tmp.atom_slice(ai,inplace=True)
        trajs.append(copy.deepcopy(tmp))
        print "TRJ "+el+" loaded"
    return trajs

def get_RMSDs(trajs,refs):
    RMSDs=[]
    for traj in trajs:
        print "new traj"
        RMSD=[]
        for ref in refs:
            print "new ref"
            traj=traj.superpose(ref,frame=0,atom_indices=ref.top.select("chainid 0"),ref_atom_indices=ref.top.select("chainid 0"))
            rmsd=calc_RMDS(traj.atom_slice(ref.top.select("chainid 1")),ref.atom_slice(ref.top.select("chainid 1")))
            RMSD.append(rmsd)
        RMSDs.append(RMSD)
    return RMSDs

def calc_RMDS(trj,ref):
    #RMSDs=[]
    coord_trj=trj.xyz
    coord_ref=ref.xyz[0]
    dist=np.square(coord_trj-coord_ref)
    RMSDs=[np.sqrt((np.sum(d))/len(d)) for d in dist]
    return RMSDs

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-top',help='topology files')
    parser.add_argument('-trajs', help='files to read', nargs='+')
    parser.add_argument('-refpdbs', help='value to consider',nargs='+')
    parser.add_argument('-out', help='output pickl file name',default='RMSD_bu.pckl', type=str)
    #parser.add_argument('-gn', help='number identifier of the group to highlight', default='428')
    return parser.parse_args()

def main():
    args = parse_args()       
    trajs = get_trjs(args.trajs,args.top) 
    refs  = get_trjs(args.refpdbs)
    RMSDs = get_RMSDs(trajs,refs)
    with open(args.out,'w') as BUF:
         pickle.dump(RMSDs,BUF)

if __name__ == '__main__': #Weird Python way to execute main()
    main()


