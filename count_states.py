#!/usr/bin/env python

import argparse
import glob as glob 
import numpy as np
import cPickle as pickle


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', help='common part of the .pckl files name',default='')
    parser.add_argument('-b', help='ignore the first b steps', default=1000,type=int)
    return parser.parse_args()

def main():
    args = parse_args()
    pickl_files=glob.glob(args.name+'*.pckl')
    res=[]
    for f in pickl_files:
        data=((np.array(pickle.load(open(f,'rb')))).transpose())[args.b:]
        Na=float(len([True for d in data if d[0]<d[1]]))/len(data)
        res.append([Na,1-Na])
    #print res
    res=(np.array(res)).transpose()
    av=np.mean(res,axis=1)
    se=np.std(res,axis=1)/np.sqrt(float(len(res)))
    for i,r in enumerate(zip(av,se)):
        print 'P({:2d}) = {:4.2f} +- {:4.2f}'.format(i,r[0],r[1]) 




if __name__ == '__main__': #Weird Python way to execute main()
    main()

