#!/usr/bin/env python

import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', help='X Y file')
    return parser.parse_args()
 
def read_XY(file_name):
    data=[]
    with open(file_name,'r') as F:
         for line in F:
             if line[0] != '@' and line[0] != '#':
                t=line.split()
                if len(t)==3:
                   print t 
                   data.append(t)
    data=np.array(data,dtype=float)     
    return(data[...,0], data[...,1])

def plot(x,y,f):
    plt.plot(x,y)
    x_new=np.linspace(x[0],x[len(x)-1],100)
    plt.plot(x_new,f(x_new))
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.show()
    return x_new, f(x_new)

def save(X,Y,file_name):
    with open(file_name+'_running_average','w') as f:
         for x,y in zip(X,Y):
            print >>f, '{:8.5f}   {:8.5f}'.format(x,y)


def main():
    args = parse_args()
    X,Y = read_XY(args.file) 
    f_fit=interp1d(X, Y, kind='cubic',bounds_error=False)
    x,y=plot(X,Y,f_fit)
    print "integral= "+str(np.trapz(y,x=x))
    #save(X,Y,args.file)

 
if __name__ == '__main__': #Weird Python way to execute main()
    main()


