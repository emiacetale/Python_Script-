#!/usr/bin/env python

import argparse
import numpy as np

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-name_in', help='Input file with E (Plum.out)',default='Plum.out')
    parser.add_argument('-name_out', help='Output file with coord (X_Plum.out)', default='X_Plum.out')
    parser.add_argument('-K', help='Spring constant (1000) [kJ/(mol nm^2)]', default=1000.0, type=float)
    parser.add_argument('-X0', help='Position of the min (0.0) [nm]', default=0.0, type=float)
    return parser.parse_args()


def getE(file_name):
    data=[]
    with open(file_name,'r') as F:
         for line in F: 
             if line[0]!='#':
                 line=line.split()
                 data.append([float(line[0]),float(line[1])])
    data=(np.array(data)).transpose()
    return(data[0],data[1])

def getXfromE(E,K,X0):
    E=np.sqrt(E*2/K)
    X=[]
    for el in E:
        if np.random.random()>0.5:
           X.append(X0+el)
        else:
           X.append(X0-el)
    return np.array(X)

def write_out(X,Y,name_out):
    with open(name_out,'w') as out:
         for x,y in zip(X,Y):
           line='{:5.3f}   {:5.3f}\n'.format(x,y)
           out.write(line)

def main():
    args = parse_args()
    T,E=getE(args.name_in)
    X=getXfromE(E,args.K,args.X0)
    write_out(T,X,args.name_out)

if __name__ == '__main__': #Weird Python way to execute main()
    main()

