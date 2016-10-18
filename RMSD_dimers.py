#!/usr/bin/env python  
 
import argparse
import mdtraj as md 


# f.write('{:5s}   {:8.6f}\n'.format(el[0],el[1]))

def get_trjs(trj_files,top='XXX'):
    trajs = []
    for el in trj_files:
        if top != 'XXX':
           trajs.append(md.load(el,top=top))
        else:
           trajs.append(md.load(el))
    return trajs




def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-top',help='topology files', nargs=1)
    parser.add_argument('-trajs', help='files to read', nargs='+')
    parser.add_argument('-refpdbs', help='value to consider',nargs='+')
    #parser.add_argument('-gn', help='number identifier of the group to highlight', default='428')
    return parser.parse_args()

def main():
    args = parse_args()       
    trajs = get_trjs(args.trajs,args.top) 
    refs  = get_trjs(args.refpdbs)
    print trajs
    print refs

    #result=files_read_and_sort(args.f,args.val) #result=[{gr1:{q1:val1_g1, q2:val2_g1, ... }, gr2:{...}}_file1, {}_file_2,...]
    #av_results=average(result)
    #plot_results(av_results, args.gn)
    #print_results(av_results,'ranking')

if __name__ == '__main__': #Weird Python way to execute main()
    main()


