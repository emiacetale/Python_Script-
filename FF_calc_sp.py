#!/usr/bin/env python

#import mdtraj as md
import argparse
import numpy as np
import subprocess as sp
import matplotlib as mpl    #To create graph directly in the queue 
mpl.use('Agg')              #"   "     "     "        "  "   " 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time 

top_num ='{}_{:02d}_{:02d}.top'

def Eone(dir, liq_name, sea_name, vac_name):
    E_LIQ=np.loadtxt(dir+'/'+liq_name,dtype=float)
    E_SEA=np.loadtxt(dir+'/'+sea_name,dtype=float)
    E_VAC=np.loadtxt(dir+'/'+vac_name,dtype=float)
    return E_LIQ/4.18400, E_SEA/4.18400 , E_VAC/4.18400 #We convert from kJ to kcal !!!

def dirs_names(charges,sigmas,F):
    l=[]
    for q in charges:
        for s in sigmas:
            l.append(F.format(s,q))
    return l

def Eall(dirs, liq_name, sea_name, vac_name):
    E_l=[]
    E_v=[]
    for dir in dirs:
        el, es, ev = Eone(dir, liq_name, sea_name, vac_name)
        E_l.append(el+es)
        E_v.append(el)
    return E_l, E_v

def DG0(E, nbin):
    P, bins = np.histogram(E, nbin, normed=1)
    bins=bins[:-1]
    dE=(bins[1]-bins[0])
    DG=np.sum(((bins+0.5*dE)*P)*dE)
    return DG    

def DDG(Ea,E0,kBT, nbin):
    ddg=[]
    for E in Ea:
        P, bins = np.histogram(E-E0,nbin)
        de=bins[1]-bins[0]
        e=-kBT*np.log(np.sum(de*P*np.exp(-(bins[:-1]+0.5*de)/kBT)))
        ddg.append(e)
        #print e 
    return ddg

def pint_res(dg_all,qrange,srange, dgt):
    DG=dg_all.reshape(len(qrange),len(srange)) # { {DG_s0, DG_s1, ... }@q0 , {DG_s0, DG_s1, ... }@q1 , ... }
    with open("results",'w') as out:
      print >>out, "#q_scale  s_scale  DG*   DG*-target"
      for q, el in zip(qrange, DG):
          for s, dg in zip(srange, el):
	      print >>out ,'{:4.2f} {:4.2f} {:9.6f} {:9.6f}'.format( q, s, dg, dg-dgt)   

def plot_res(dg, qrange, srange):
    m=np.amin(dg)
    M=np.amax(dg)
    if (-m)>M: M=-m
    dg=dg.reshape(len(qrange),len(srange)) 
    dg=np.transpose(dg)
    dq=qrange[1]-qrange[0]
    qrange=qrange-0.5*dq
    qrange=np.append(qrange, qrange[-1]+dq)
    ds=srange[1]-srange[0]
    srange=srange-0.5*ds
    srange=np.append(srange, srange[-1]+ds)
    fig=plt.figure()
    ax1 = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    ax10= fig.add_axes([0.91, 0.15, 0.02, 0.75])

    #print srange
    #print qrange
    #print dg
    
    ax1.pcolor(qrange, srange, dg, vmin=-M, vmax=M, cmap='bwr')
    ax1.set_xlabel(r'$\Delta q$',fontsize=25)
    ax1.set_ylabel(r'$\Delta \sigma$',fontsize=25)
    ax1.set_xlim(srange[0],srange[-1])
    ax1.set_ylim(qrange[0],qrange[-1])
    ax1.set_xticks(np.arange(qrange[0]+0.5*dq,qrange[-1],dq*5))
    ax1.set_yticks(np.arange(srange[0]+0.5*ds,srange[-1],ds*5))
    ax1.tick_params(axis='y', which='major', labelsize=20)
    ax1.tick_params(axis='x', which='major', labelsize=20)

    cb=mpl.colorbar.ColorbarBase(ax10, cmap='bwr', orientation='vertical', norm=mpl.colors.Normalize(vmin=-M, vmax=M),ticks=[-M, -M/2, 0, M/2, M])
    plt.savefig('Result.png')

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-dgt',type=float, help='target solvation free energy (kcal)')
    parser.add_argument('-np', type=int, help='number of points [21]', default=21)
    parser.add_argument('-qs_min', type=float, help='min scaling of q [0.8]',default=0.8)
    parser.add_argument('-qs_max', type=float, help='max scaling of q [1.2]',default=1.2)
    parser.add_argument('-ss_min', type=float, help='min scaling of s [0.8]', default=0.8)
    parser.add_argument('-ss_max', type=float, help='max scaling of s [1.2]', default=1.2)
    parser.add_argument('-nbin',type=int, help='number of bins to calulate P(E), P(DE) distr [50]', default=50)
    parser.add_argument('-e_liq_name', type=str, help='file of liq intra-E [E_LIQ]', default='E_LIQ')
    parser.add_argument('-e_sea_name', type=str, help='file of liq inter-E (SEA) [E_SEA]', default='E_SEA')
    parser.add_argument('-e_vac_name', type=str, help='file of vac intra-E [E_VAC]', default='E_VAC')
    parser.add_argument('-format_dir', type=str, help='foramt dir name [{:4.2f}]', default='{:4.2f}')
    return parser.parse_args()

def main():
    #"Hard-coded" variables: not much sense to change them at this point, but it is nice to have an handle for future 
    kBT=300*1.9872041E-3   # KT in kcal/mol

    #Script:
    args = parse_args()                                                                      #We read some input 
    E0_l, E0_s, E0_v=Eone((args.format_dir+'_'+args.format_dir).format(1.0,1.0), 
         args.e_liq_name, args.e_sea_name, args.e_vac_name) 
    dirs=dirs_names(np.linspace(args.qs_min, args.qs_max, num=args.np),
         np.linspace(args.ss_min, args.ss_max, num=args.np),args.format_dir+'_'+args.format_dir)
    E_l, E_v = Eall(dirs,args.e_liq_name, args.e_sea_name, args.e_vac_name)
    dg0=DG0(E0_s,args.nbin)                               # DG* original FF
    ddg_L= np.array(DDG(E_l,E0_l+E0_s,kBT,args.nbin))     # DDG all liq

    ddg_V= np.array(DDG(E_v,E0_v,kBT,args.nbin))          # DG all vac
    dg_all=(ddg_L-ddg_V)+dg0                              # DG*=DDG+DG0*-DGv
    pint_res(dg_all, np.linspace(args.qs_min, args.qs_max, num=args.np), 
         np.linspace(args.ss_min, args.ss_max, num=args.np), args.dgt)          # Print results
    plot_res(dg_all-args.dgt, np.linspace(args.qs_min, args.qs_max, num=args.np), 
         np.linspace(args.ss_min, args.ss_max, num=args.np))           # Plot results

if __name__ == '__main__': #Weird Python way to execute main()
    main()



