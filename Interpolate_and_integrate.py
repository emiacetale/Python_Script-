#!/usr/bin/env python

import numpy as np
from scipy import interpolate
import argparse
import matplotlib.pyplot as plt


def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='input x y dy file')
    return parser.parse_args()

def main():
    args = parse_args()
    x,y,dy=np.loadtxt(args.f, unpack=True)
    f = interpolate.interp1d(x,y,kind='cubic',assume_sorted=False, bounds_error=False)
    xnew=np.arange(0.0,1.01,0.01)
    ynew= f(xnew)
    #plt.plot(x, y, 'o', xnew, ynew, '-')
    #plt.show() 
    #print xnew
    #print ynew
    Int_n=np.trapz(ynew,xnew)
    Int=np.trapz(y,x)
    print Int_n, Int
 
if __name__ == '__main__': #Weird Python way to execute main()
    main()

