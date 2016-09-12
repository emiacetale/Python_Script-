#!/usr/bin/env python

import argparse
#import glob as glob 
import numpy as np
import cPickle as pickle

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', help='common part of the files name [dhdl]',default='dhdl')
    parser.add_argument('-nlambda', help='number of lambda points [ 101==[0:100] ]', default=101, type=int)
    parser.add_argument('-lambda_step', help='consider only every lambda step (for test purpose)', default=1, type=int)
    parser.add_argument('-nsim', help='number of sim for lambda point [ 50==[0:49] ]', default=50, type=int)
    parser.add_argument('-b', help='ignore the first b steps', default=1000,type=int)
    parser.add_argument('-graph',help='plot the graph of dhdl vs lambda', default=False, action='store_true')
    parser.add_argument('-ba',help='calc error using the block av of n blocks [0 == std error]', default=0, type=int)
    return parser.parse_args()

def main():
    args = parse_args()
    dhdl=[]
    lambda_step=1.0/float(args.nlambda-1)
    for i in np.linspace(0, args.nlambda, num=args.nlambda/args.lambda_step, endpoint=False, dtype=int):
        #l_t,dhdl_t=read_lambda(args.)
        #print i
        dhdl.append(av_dhdl(args.name+"_"+str(i),i,args.nsim,args.b, args.ba))
    dhdl=np.array(dhdl).transpose()
    with open('dhdl_data.pckl','wb') as Fout:
         pickle.dump(dhdl,Fout)
    #print dhdl
    DG=np.trapz(dhdl[0],dx=lambda_step)
    error=calc_error_integral(dhdl[1],lambda_step)
    print "DG = ",DG," +- ", error," kJ/mol"
    if args.graph:
        plot_graph(dhdl,args.nlambda)

def calc_error_integral(data,h):
    error=[ 0.5*h*np.sqrt(a*a+b*b) for a,b in zip(data[1:],data[:-1]) ] #Error of each trapezoid
    return np.sqrt(np.sum(np.square(error)))                            #"Sum" of the trapezoid errors

def av_dhdl(name_head,lambda_point,nsim,step_begin,blocks):
    lambda_point=lambda_point+1 #We are going to read the column lp+1, since the first one is the step
    dhdl=[]
    for i in np.linspace(0, nsim, num=nsim, endpoint=False, dtype=int):
        #dhdl.append(load_data(name_head+'_'+str(i)+'.xvg',lambda_point,np.array(['@','#'])))
        dhdl.append(load_data(name_head+'_'+str(i)+'.xvg',1,np.array(['@','#'])))
    dhdl=(np.array(dhdl).flatten())[step_begin:]
    if blocks==0:
        error=(np.std(dhdl)/np.sqrt(len(dhdl)))
    else:
        error=block_av(dhdl,blocks)
    return([np.average(dhdl),error])

def block_av(data,n_blocks):
    #print len(data), n_blocks, len(data)/n_blocks
    data=np.array_split(data,n_blocks)
    block=np.array([np.average(d) for d in data])
    return np.std(block)/np.sqrt(n_blocks-1)

def load_data(filename, column, comments_characters):
    out=[]
    with open(filename,'r') as F:
         for line in F: 
             if not (line[0]==comments_characters).any(): #Checks that the line is not a comment line
                out.append(float(line.split()[column]))
    return out

def plot_graph(dhdl,lambda_step):
    import matplotlib as mpl    #To create graph directly in the queue 
    mpl.use('Agg')              #"   "     "     "        "  "   "   
    import matplotlib.pyplot as plt 
    #import matplotlib.cm as cm
    import seaborn as sns 
    
    #y=dhdl[0]
    x=np.linspace(0, 1, lambda_step)
    sns.set(style="darkgrid", palette="Set2")
    plt.plot(x,dhdl[0])
    plt.errorbar(x,dhdl[0],yerr=dhdl[1],fmt='')
    #plt.show()
    plt.savefig('dHdl_vs_lambda.png')

if __name__ == '__main__': #Weird Python way to execute main()
    main()

