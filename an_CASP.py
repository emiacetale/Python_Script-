#!/usr/bin/env python  
 
import argparse
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from copy import deepcopy


def file_operations(f,val_list):                   #input: file to read , value we want to bring back
    results=[]
    with open(f) as input:
         lines=input.readlines()
    lines[0]=lines[0].split()                      #will be the key of the dict
    for line in lines[1:]:
      if len(line)>1:
        rg={}
        line=line.split()
        #print line
        line[1]=(line[1].split('_'))[1]            #We keep only the model submission number
        for name,el in zip(lines[0],line):
            rg.update({name:el})
        results.append(rg)                         #Array of dicts with all the results
    #print results
    interesting={}
    for el in results:                             #We now create a dict of dicts {Group_num:{value1:xxx,value2:yyy}}
        interesting.update({el['GR#']:{key:el[key] for key in val_list}})
    return interesting

def files_read_and_sort(files,val_list):  #list of files, list of value we are interested in 
    container=[]
    for f in files: 
        container.append(file_operations(f,val_list))
                                          #Now container is an array of dictionary whom keys are the group number
                                          #each has a dict with all the quantityes we are interested in 
    container=clean_result(container)     #keeps only the groups that predicted ALL the structures (fair comparison???)
    #print container[0].keys()
    return container

def clean_result(container):              #container=[{G1:{quantity1:value1},{G2:{...},...}}_struct1,{}_struct2,...]
    #print array[0]
    data=set(container[0].keys())         #we look for the groups that are in every structure  
    for el in container[1:]:              #   using set and intersection between sets 
        data=(data & set(el.keys()))      #
    new_container=[]                                  
    for el in container:
        new_container.append({key:el[key] for key in data})  #We keep only the groups that predicted all the structures we passed
    return new_container

#def plot_results(result_container):       #we plot the average of each group for each quantity considered 
 
def average(result_container):
    av_container=deepcopy(result_container[0])#will contain the average value of every quantity for every group
    n_entry=deepcopy(result_container[0]) #somtimes we have N/A so we need to average on a different number of entries 
    group_keys=n_entry.keys()             #All the groupnames 
    quantity_keys=n_entry[group_keys[0]].keys()  #all the quantity names
    #print group_keys, quantity_keys
    for group in group_keys:              #We initialize the 2 collector of the results 
        for quantity in quantity_keys:
            n_entry[group][quantity]=0
            av_container[group][quantity]=float(0.0)
    #print av_container
    #print n_entry
    for struct in result_container:       #We collect the values...  
        for group in group_keys:
            for quantity in quantity_keys:
                if struct[group][quantity] != 'N/A':  #avoid problem with N/A
                   #print struct[group][quantity]
                   n_entry[group][quantity]=n_entry[group][quantity]+int(1)
                   av_container[group][quantity]=av_container[group][quantity]+float(struct[group][quantity])
                   #print n_entry[group][quantity], av_container[group][quantity]
    for group in group_keys:              #Calculates the averages          
        for quantity in quantity_keys:
            av_container[group][quantity]=av_container[group][quantity]/float(n_entry[group][quantity])
    #print av_container
    return av_container

def plot_results(results, group_name):
    group_keys=results.keys()
    quantity_keys=results[group_keys[0]].keys()
    for quantity in quantity_keys:
        data=[]
        for group in group_keys:
            data.append(float(results[group][quantity]))
        data=np.sort(np.array(data))
        data_us=float(results[group_name][quantity]) 
        plot_result(data,data_us,quantity,group_name)

def plot_result(data,data_us,quantity_name,group_name):
    bar_size=0.2
    fig=plt.figure()
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    all_bar = ax.bar((np.arange(data.size))*bar_size, data, width=0.2,color='gray')
    us_bar = ax.bar(0.2*(np.where(data==data_us))[0], [data_us], width=0.2,color='orange')
    ax.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax.tick_params(axis='y', labelsize=20)
    ax.set_ylabel(quantity_name,fontsize=25)
    ax.set_xlim(0.0,data.size*bar_size)
    plt.savefig('Result_'+quantity_name+'.png')

def print_results(results,name_string):
    group_keys=results.keys()
    quantity_keys=results[group_keys[0]].keys()
    for quantity in quantity_keys:
        data=np.array([[group_keys[0],float(results[group_keys[0]][quantity])]],dtype=object)
        for group in group_keys[1:]:
            data=np.append(data,[np.array([group,float(results[group][quantity])], dtype=object)],axis=0)
        data=np.sort(data,axis=0)
        print_result(data,name_string+'_'+quantity)

def print_result(data,filename):
    with open(filename,'w') as f:
         for el in data:
             f.write('{:5s}   {:8.6f}\n'.format(el[0],el[1]))

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    default_string=['GDT_TS','GDT_HA','RMS_CA','RMS_ALL','RMSD[L]','RMSD[D]','MolPrb_Score','SphGr','FlexE','TMscore']
    parser.add_argument('-f', help='files to read', nargs='+')
    parser.add_argument('-val', help='value to consider',nargs='+',default=default_string)
    parser.add_argument('-gn', help='number identifier of the group to highlight', default='428')
    return parser.parse_args()

def main():
    args = parse_args()        
    result=files_read_and_sort(args.f,args.val) #result=[{gr1:{q1:val1_g1, q2:val2_g1, ... }, gr2:{...}}_file1, {}_file_2,...]
    av_results=average(result)
    plot_results(av_results, args.gn)
    print_results(av_results,'ranking')

if __name__ == '__main__': #Weird Python way to execute main()
    main()


