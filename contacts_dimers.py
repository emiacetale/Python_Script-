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
        ai=tmp.top.select("name CB")
        tmp.atom_slice(ai,inplace=True)
        trajs.append(copy.deepcopy(tmp))
        print "TRJ "+el+" loaded"
    return trajs

def get_pairs(struct):
    list1=struct.top.select("chainid 0")
    list2=struct.top.select("chainid 1")
    pairs=[ [i,j] for i in list1 for j in list2 ]
    #print pairs
    return np.array(pairs)

def get_contacts(trajs,pairs_list,dist):
    CONTACTS=[]
    for traj in trajs:
        contact=get_contact(traj,pairs_list,dist)
        CONTACTS.append(contact)
    return CONTACTS

def get_contact(traj,pairs_list,dist):
    contact=(md.compute_distances(traj, pairs_list, periodic=False))<dist
    return contact


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-top',help='topology files')
    parser.add_argument('-trajs', help='files to read', nargs='+')
    parser.add_argument('-contact_dist', help='Atoms closer than this distance [nm] are in contact [0.8]',type=float,default=0.8)
    parser.add_argument('-refpdbs', help='value to consider',nargs='+')
    parser.add_argument('-out_traj', help='prefix for traj results',default='contacts_', type=str)
    parser.add_argument('-out_ref', help='output for ref results pickl file name',default='Contacts_ref_bu.pckl', type=str)
    return parser.parse_args()

def main():
    args = parse_args()       
    trajs = get_trjs(args.trajs,args.top) 
    refs  = get_trjs(args.refpdbs)
    pairs_list=get_pairs(refs[0]) #we get the pair list from the reference structrure
    contacts_refs=get_contacts(refs,pairs_list,args.contact_dist) #save contacts of refs
    with open(args.out_ref,'w') as BUF:
             pickle.dump(contacts_refs,BUF)
    for traj,name in zip(trajs,args.trajs):
        contacts=get_contact(traj,pairs_list,args.contact_dist)
        with open(args.out_traj+name+'.pckl','w') as BUF:
             pickle.dump(contacts,BUF)


             

if __name__ == '__main__': #Weird Python way to execute main()
    main()


