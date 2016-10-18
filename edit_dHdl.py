#!/usr/bin/env python

import glob as glob
import numpy as np
import argparse
import matplotlib.pyplot as plt


def get_file_list_and_el(pref,suf,skip):  
    l=[[e[0],int(e[1]),e[2],e[3]] for e in [el.split('.') for el in glob.glob(pref+"*"+suf)]]
    l.sort(key=lambda el: el[1])
    select=np.linspace(0,len(l)-1,num=int(len(l)/skip),dtype=int)
    return([[el[0]+"."+str(el[1])+"."+el[2]+"."+el[3] for el in l][i] for i in select], select)

def readfile(filename):
    h1=[]
    h2=[]
    d=[]
    with open(filename,'r') as F:
         for line in F:
            if line[0]=='#': 
                h1.append(line)
            elif line[0]=='@': 
                h2.append(line)
            else:
             d.append(line)
    return (h1, h2, d)     
         
def edit_h2(lines,sel):
    h=[]
    b=[]
    for l in lines:
        #print l
        if l[2]=='s' and l[6]=='l':
            b.append(l)
            print l 
        else:
            h.append(l)
    b=[b[i] for i in sel]
    for l in b:
        h.append(l)
    return(h)

def edit_data(data,sel):
    d=np.array([el.split() for el in data])
    d=d.transpose()
    d=[d[i] for i in sel]
    d=d.transpose()
    out=[]
    for line in d:
        s=''
        for el in line:
            s=s+(str(el))+" "
        out.append(s)
    return(out)

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', help='file prefix',default='out')
    parser.add_argument('-q', help='file suffix',default='dhdl.xvg')
    parser.add_argument('-skip', help='read only every nfiles', default=2, type=int)
    return parser.parse_args()
 
def main():
    args = parse_args()
    files_list, sel=get_file_list_and_el(args.p,args.q,args.skip)
    for f in files_list:
        head1, head2, data = readfile(f)
        #print head2
        head2=edit_h2(head2,sel)
        #data=edit_data(data,sel)

 
if __name__ == '__main__': #Weird Python way to execute main()
    main()


