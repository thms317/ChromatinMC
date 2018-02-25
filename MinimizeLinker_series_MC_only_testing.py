# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:34:39 2016

@author: Visscher
"""

import FiberMC as FMC
import NucleosomeMC as NMC
import numpy as np
from helixmc import random, kBT
from helixmc.pose import HelixPose
import time


def GetCut(DNA, nucl, n_nuc, cut, dyads):
     fiber_params = FMC.set_fiber_params(n_nuc, diameter=330.0, rise=100.0, nld=17, ribs=1)
     # fiber_params[0,5] = 1.6
     print(fiber_params)
     fiber_dyad_frames = FMC.get_fiber_dyad_frames(fiber_params, nucl.d_index, nucl.nuc_type)
     print (fiber_dyad_frames)     
     FMC.cast_nucs_fiber(DNA, dyads, fiber_dyad_frames, cut=cut)




def MC(DNA, start=179, end=239, unwrap=-20):
    for i in range(start + unwrap, end - unwrap):

        Estart = FMC.Boltzmann_Energy(DNA.params[start + unwrap]) + FMC.Boltzmann_Energy(DNA.params[i])
        pars = DNA.params
        pars[i] += (0.5 - np.random.rand()) * np.array([0.001, 0.001, 0.01, 0.001, 0.001, 0.01])
        newDNA = DNA.copy()
        newDNA.set_params(pars)
        GetCut(newDNA, nucl, 2, start - dyads[0] + unwrap, dyads)
        Enew = FMC.Boltzmann_Energy(newDNA.params[start + unwrap]) + FMC.Boltzmann_Energy(newDNA.params[i])

        if Enew - Estart < 0:
            DNA.set_params(newDNA.params)
            print Estart - Enew


def FindMinimum(DNA, nucl, a, dyads, start=177, end=239, unwrap=-20):
    Eold = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
    minimum = True
    for i in range(5):
        print i
        Eprev = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
        a[i] += 0.001
        CurveDNA(a, unwrap, dna, nucl, start, end, 2, dyads)
        Ecurr = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
        while Ecurr < Eprev:
            minimum = False
            Eprev = Ecurr
            a[i] += 0.001
            CurveDNA(a, unwrap, dna, nucl, start, end, 2, dyads)
            Ecurr = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
            print '+'
        if minimum == False:
            Eprev = Ecurr
            a[i] -= 0.001

        if minimum:
            a[i] -= 0.002

        CurveDNA(a, unwrap, dna, nucl, start, end, 2, dyads)
        Ecurr = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
        while Ecurr < Eprev:
            minimum = False
            Eprev = Ecurr
            a[i] -= 0.001
            CurveDNA(a, unwrap, dna, nucl, start, end, 2, dyads)
            Ecurr = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
            print '-'
        a[i] += 0.001
        CurveDNA(a, unwrap, dna, nucl, start, end, 2, dyads)

        print Eprev, Ecurr

    Enew = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
    print a
    print Enew, Eold
    return a, Ecurr, minimum


def CurveDNA(apars, unwrap, DNA, nucl, start, end, n_nuc, dyads):
    """
    Curve the linker DNA
    """
    newpars = DNA.params

    newpars[start + unwrap:end - unwrap, :] = FMC.create_curved_linker(end - start - 2 * unwrap, apars)
    DNA.set_params(newpars)
    GetCut(DNA, nucl, n_nuc, start - dyads[0] + unwrap, dyads)


# ititiate program

min_unwrap = 0
max_unwrap = 1

date = time.strftime("%Y%m%d")
NRL = 197
file_name = 'simulations\\NLD test\\' + date + '_' + str(NRL)

read_from_file = False

n_nuc = 2
n_bp = 147 + NRL
dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL)
# end = FMC.nuc_get_ends(dna, dyads[1], nucl)[0]
start = dyads[0] + (147 - nucl.d_index)
end = dyads[1] - nucl.d_index
# start = FMC.nuc_get_ends(dna, dyads[0], nucl)[1]

if read_from_file == True:

    dna_file = 'simulations\\20180215_197_0_unwrap.npz'
    dna = HelixPose(params=np.load(dna_file)['params'], frame0=np.load(dna_file)['frame0'])

else:

    # 197
    a = np.array([11.617, 0.329, 3.672, 0.965, 0.800])  # 15/02/18 197 solenoid --> unwrap 0-10
    # a = np.loadtxt('simulations\\180215_197_0_unwrap.dat')

    unwrap = 0

    link_len = end - start

    CurveDNA(a, unwrap, dna, nucl, start, end, n_nuc, dyads)
    print start
    print end
    se = np.linalg.norm(dna.coord[start] - dna.coord[end])



    # def GetDot(DNA, start, end):
    #    A1 = DNA.frames[start+unwrap]
    #    A2 = DNA.frames[end-unwrap]
    #    z = np.array([0,0,1])
    #    return np.dot(np.dot(A1,z), np.dot(A2,z))
    #
    # theta = np.arccos( GetDot(dna, start, end) )
    # print theta*se/(2*np.sin(.5*theta))
    #
    Ecurr = FMC.Boltzmann_Energy(dna.params[start + unwrap])
    
    '''
    
    
    print dna.params[start + unwrap]
    print Ecurr
    minimum = False

    # optimize the Bert/John parameters
    while minimum == False:
        a, Ecurr, minimum = FindMinimum(dna, nucl, a, dyads, start=start, end=end, unwrap=unwrap)

    Etemp = FMC.Boltzmann_Energy(dna.params[start + unwrap])
    print Etemp

    dna.write2disk(str(file_name) + '_' + str(unwrap) + '_unwrap_optimized')
    np.savetxt(str(file_name) + '_' + str(unwrap) + '_unwrap_input_params_optimized.dat', a)

E_res = np.array([])


for unwrap in range(min_unwrap,max_unwrap+1):

    unwrap = -unwrap

    E = []
    j = 0

    print start
    print end
    number_of_iterations = 3
    while j < number_of_iterations:
        MC(dna, start, end, unwrap)
        print str(j+1) + ' of ' + str(number_of_iterations) + ' (' + str(unwrap) + ' base pairs unwrapped)'
           # nice to include the energy calculated from the helixmc score-function
        Etemp = FMC.Boltzmann_Energy(dna.params[start + unwrap])
        print Etemp
        E.append(Etemp)
        j += 1

    np.savetxt(str(file_name) + '_' + str(unwrap) + '_unwrap_E.dat ', E)
    dna.write2disk(str(file_name) + '_' +str(unwrap)+'_unwrap')

    E_res = np.concatenate((E_res, E))

E_res = np.transpose(np.hsplit(E_res,max_unwrap+1))
np.savetxt(str(file_name) + '_energies_series.dat',E_res)

    '''