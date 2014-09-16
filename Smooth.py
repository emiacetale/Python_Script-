#!/usr/bin/env python

import numpy as np
import argparse


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', help='X Y file')
    parser.add_argument('-n',    help='number of time to smooth [10]',type=int, default=10)
    return parser.parse_args()
 
def read_XY(file_name):
    data=np.loadtxt(file_name,dtype=float)
    return(data)

def smooth(data,n):
    dx=(data[1,0]-data[0,0])
    data=np.insert(data, 0, [[(data[0,0]-dx),0]], axis=0)
    data=np.append(data,[[(data[len(data)-1,0]+dx),0]],axis=0)
    for interation in range(0,n):                                        
       data[1:-1,1]=0.5*data[1:-1,1]+0.25*data[0:-2,1]+0.25*data[2:,1] #This smart move is called slicing (better than loop) 
    data=data[1:-1]
    data[:,1]=data[:,1]/np.trapz(data[:,1], data[:,0])
    return(data)    

def print_value_raw(value,Fname,formatting):            #prints a matrix of values (lenght of lines can be variable)
    formatting=formatting+"   "
    with open(Fname,'w') as F:
        for el in value: 
            for n in el: 
                F.write(formatting.format(n))
            F.write("\n")
        

def main():
    args = parse_args()
    data = read_XY(args.file) 
    data_new = smooth(data,args.n)
    print_value_raw(data_new,args.file+'_NEW','{:8.6f}')


 
if __name__ == '__main__': #Weird Python way to execute main()
    main()


