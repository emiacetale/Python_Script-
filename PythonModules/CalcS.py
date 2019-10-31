#This is the module that is being sourced to compute Information 

import numpy as np 
from scipy.special import digamma as dg

def Calc_s(n,N):
    if n==0:
        return 0 
    else:
        return (n/N)*(np.log(N)-dg(n)-(np.power(-1,n)/(n+1.0))) 
    
def Calc_S(data):
    N=float(len(data))
    data=np.array(data,dtype=float)
    #print N, np.max(data)/N, np.power(-1,np.max(data)), np.max(data)
    return np.sum([Calc_s(ni,N) for ni in data ])
        
def CalcI(px,py,pxy):
    return Calc_S(px)+Calc_S(py)-Calc_S(pxy.flatten())

class Info:
    def __init__(self, x, y, nbin=24):
        self.nb=nbin #We want to have it saved somewhere 
        self.px, edg  = np.histogram(x,bins=nbin,density=False)
        self.val_px= [e+0.5*(edg[0]+edg[1]) for e in edg[:-1]]
        self.py, edg  = np.histogram(y,bins=nbin,density=False)
        self.val_py= [e+0.5*(edg[0]+edg[1]) for e in edg[:-1]]
        self.pxy, xed, yed = np.histogram2d(x,y,bins=nbin,density=False)
        #We don't need to save xed and yed since they should be the same of 1D hystograms 
        #self.val_pxy= [ [e+0.5*(xed[0]+xed[1]) for e in xed[:-1]] , [e+0.5*(yed[0]+yed[1]) for e in yed[:-1]] ]
        self.I=CalcI(self.px,self.py,self.pxy)



