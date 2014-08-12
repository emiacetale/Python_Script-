#! /usr/bin/env python

import sys
import numpy
from pyevolve import G1DList, GSimpleGA, Initializators, Mutators, Consts, Selectors
from pyevolve import Scaling

# Uses Use the Lifson - Roig Model to find helical propensities of different amino acids.
# J. Phys Chem B,2009,113,9004-9015
# I (Emiliano) am trying to optimize a bit Alberto scripts in a single one 
# Round 1 refers to the analysis of pure alanine pept
# Round 2 refers to the analysis of Best's pept (AAXAAAAXAAAAXAA) 

n_steps = 100                    #This has to go inside main!       
#v_ALA = 0.54                     # "   "   "  "  "      "
#w_ALA = 1.36                     # "   "   "  "  "      "

#After run 10
#v_ALA = 0.2528024035614776       #These are the output of the 1st round and the input of the second (this also has to go in main)
#w_ALA = 1.3932075952692177       #This has to go inside main! 

def tovw(mat):                   #Same in both Al's scripts  
    '''
    make a logical matrix:  
        0 --> coil
        1 --> helical segment
        -1--> helical
    '''
    my_s = mat.shape
    h_mat = numpy.zeros( my_s )
    for i in range(len(mat)):
        if mat[i,0] == 1:
            h_mat[i,0] = -1
        if mat[i,len(mat[i])-1] == 1:
            h_mat[i,len(mat[i])-1] = -1
        for j in range(1,len(mat[i])-1):
            if mat[i,j] == 0:
                continue
            if (mat[i,j-1] == 0) or (mat[i,j+1] == 0):
                h_mat[i,j] = -1
                continue
            h_mat[i,j] = 1
    return h_mat

def make_partition1(res,v=10,w=100):                      #Same in both
    r1 = [w,v,0]
    r2 = [0,0,1]
    r3 = [v,v,1]
    m  = [r1,r2,r3]
    m = numpy.matrix(m)
    c = numpy.matrix(numpy.identity( 3 ))
    for i in range(res):
        c = c * m
    a =  numpy.matrix([0,1,1])
    b = numpy.matrix([0,0,1])
    c = b* c * a.T
    return numpy.matrix.item(c)

#def make_Best_partition(res,v=10,w=100):
def make_partition(res,seq,w_ALA,v_ALA,v=10,w=100):   #(res, 'AAAAAAAAAAAAAAA', 1, 1, 10, 100) for first use 
                                                    #(res, 10, 100, 'AAXAAAAXAAAAXAA', w_ALA, v_ALA) for second use   

    r1 = [w,v,0]                            # This creates the initial condition for both round
    r2 = [0,0,1]                            #                    
    r3 = [v,v,1]                            #
    m  = [r1,r2,r3]                         #
    m = numpy.matrix(m)                     #

    r1_ALA = [w_ALA,v_ALA,0]                # This is important only on the second round, since
    r2_ALA = [0,0,1]                        # it will not contributre to the matrix c until there
    r3_ALA = [v_ALA,v_ALA,1]                # is only ALA in the sequence. 
    m_ALA  = [r1_ALA,r2_ALA,r3_ALA]         # Nevertheless we need to create this both time 
    m_ALA = numpy.matrix(m_ALA)             #

    c = numpy.matrix(numpy.identity( 3 ))   # The first round (only ALA pept) the seq will be 
    for Aa in seq:                          # 'AAAAAAAAAAAAXAA' therefore c is equal to the c 
        if "A" in Aa:                       # created in Alberto's script LR.py
            c = c * m_ALA                   #
        if "X" in Aa:                       # The second round the seq will be different,this
            c = c * m                       # gives back situation of Al's LR_mix.py 

    a =  numpy.matrix([0,1,1])
    b = numpy.matrix([0,0,1])
    c = b* c * a.T

    return numpy.matrix.item(c)

#def scoring_function_factory(Nconf,Nw,Nv,Nw_x,Nv_x,N_res):
def scoring_function_factory(Nconf,Nw,Nv,Nw_x,Nv_x,N_res,seq,w_ALA,v_ALA): #R1 (Nconf,0,0,Nw_x,Nv_x,N_res,'AA..AA',1,1)
                                                                           #R2 ()
    def eval(chromosome):
        v = chromosome[0]
        w = chromosome[1]*10
        Z = make_partition(N_res,seq,w_ALA,v_ALA,v=v,w=w)
        #Z = make_partition1(N_res,v=v,w=w)
        
        #R1:        0                             0                 num                 num                   num
        #R2:        num                           num               num                 num                   num
        lnL = Nw * numpy.log(w_ALA) + Nv * numpy.log(v_ALA) + Nw_x * numpy.log(w) + Nv_x * numpy.log(v) - Nconf * numpy.log(Z)
        '''Will be a negative value, want this number to increase,multiply by -1 and minimize instead of maximize (the algorithm cannot handle negative scores)'''
        #print lnL
        return -lnL
    return eval

def prep_A(psi,phi):
    l_phi = (phi <= -30.) * (phi > -100.) * 1
    l_psi = (psi <= -7.) * (psi >= -67.) *1
    h_matrix = l_phi*l_psi
    logical_mat = tovw(h_matrix)

    number_conf = len(logical_mat)
    number_helical_w = ((logical_mat > 0) * 1).sum()
    number_helical_v = ((logical_mat < 0) * 1).sum()
    return (number_conf, 0,    0,    number_helical_w, number_helical_v, len(logical_mat[0]))
           #Nconf,       Nw,   Nv,   Nw_x,             Nv_x,             N_res

def prep_P(psi,phi):
    l_phi = (phi <= -30.) * (phi > -100.) * 1
    l_psi = (psi <= -7.) * (psi >= -67.) *1
    h_matrix = l_phi*l_psi
    logical_mat = tovw(h_matrix)
    number_conf = len(logical_mat)
    number_helical_w = ((logical_mat > 0) * 1).sum()
    number_helical_v = ((logical_mat < 0) * 1).sum()
    number_helical_w_resx = ((logical_mat[:,2] > 0) * 1).sum() +  ((logical_mat[:,7] > 0) * 1).sum() +  ((logical_mat[:,12] > 0) * 1).sum()
    number_helical_v_resx = ((logical_mat[:,2] < 0) * 1).sum() +  ((logical_mat[:,7] < 0) * 1).sum() +  ((logical_mat[:,12] < 0) * 1).sum()
 
    number_helical_w -= number_helical_w_resx 
    number_helical_v -= number_helical_v_resx 
 
    n_res = len(logical_mat[0])
    #partition_f = make_partition(len(logical_mat[0]),v=1,w=1) (len(logical_mat[0]),v=1,w=1,'AAAAAAAAAAAAAAA', 1, 1)
    return (number_conf, number_helical_w, number_helical_v, number_helical_w_resx, number_helical_v_resx, len(logical_mat[0]))
           #Nconf,       Nw,               Nv,               Nw_x,                  Nv_x,                  N_res

def ga_magic(eval_func, n_steps):
    genome = G1DList.G1DList(2)
    genome.evaluator.set(eval_func)
    genome.setParams(rangemin=0.000001, rangemax=1.0, gauss_sigma=0.01)
    genome.initializator.set(Initializators.G1DListInitializatorReal)
    genome.mutator.set(Mutators.G1DListMutatorRealGaussian)
    ga = GSimpleGA.GSimpleGA(genome)
    ga.setMinimax( Consts.minimaxType["minimize"] )
    ga.setElitism(True)
    #ga.selector.set(Selectors.GRouletteWheel)
    ga.selector.set(Selectors.GRankSelector)
    ga.nGenerations = n_steps
    ga.setMutationRate(0.20)
    ga.evolve(freq_stats=0) #Was 10
    res=ga.bestIndividual()
    return res[0], res[1]*10

def eval_func(chromosome):
   score = 0.0
   # iterate over the chromosome
   for value in chromosome:
      if value==0:
         score += 1
   return score


def main(f_phi_A,f_psi_A,f_phi_P,f_psi_P):
 
    #AAAAAAAA part
    
    phi = numpy.loadtxt(f_phi_A)
    psi = numpy.loadtxt(f_psi_A)
 
    inp=prep_A(psi, phi)
    
    eval_func = scoring_function_factory(inp[0],inp[1],inp[2],inp[3],inp[4],inp[5],'XXXXXXXXXXXXXXX',1,1)
                                        #Because of Al's twisted logic we need to declare a pol made of all A as 
                                        #made of all X... Don't worry it is fine like that. I checked :)
    v_A,w_A=ga_magic(eval_func, 100)
    print 'ALA v and w: ',v_A, w_A

    #AAXAAAAXAAAAXAA part

    phi = numpy.loadtxt(f_phi_P)
    psi = numpy.loadtxt(f_psi_P)

    inp=prep_P(psi, phi)
    eval_func = scoring_function_factory(inp[0],inp[1],inp[2],inp[3],inp[4],inp[5],'AAXAAAAXAAAAXAA',w_A,v_A)

    v_P,w_P=ga_magic(eval_func, 100)
    print 'PEP v and w: ',v_P, w_P





if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2], sys.argv[3],sys.argv[4])

