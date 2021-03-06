#!/usr/bin/env python

import numpy as np
import argparse
import matplotlib.pyplot as plt


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', help='X Y file')
    parser.add_argument('-b',help='ps to skip [1000]',default=1000, type=float)
    return parser.parse_args()
 
def read_XY(file_name,b):
    data=[]
    with open(file_name,'r') as F:
         for line in F:
             if line[0] != '@' and line[0] != '#':
                t=line.split()
                if float(t[0])>b:
                   data.append(t)
    data=np.array(data,dtype=float)     
    return(data[...,0], data[...,1])

def running_av(data):
    s=np.cumsum(data)
    for n,el in enumerate(s):
        s[n]=el/(n+1)
    return s

def plot(x,y):
    plt.plot(x,y)
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.show()

def save(X,Y,file_name):
    with open(file_name+'_running_average','w') as f:
         for x,y in zip(X,Y):
            print >>f, '{:8.5f}   {:8.5f}'.format(x,y)


def main():
    args = parse_args()
    X,Y = read_XY(args.file,args.b) 
    Y = running_av(Y)
    plot(X,Y)
    save(X,Y,args.file)

 
if __name__ == '__main__': #Weird Python way to execute main()
    main()


