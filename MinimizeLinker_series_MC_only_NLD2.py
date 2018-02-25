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
import os


def GetCut(DNA, nucl, n_nuc, cut, dyads, diameter=330.0, rise=100.0, nld=100, ribs=1):
    fiber_params = FMC.set_fiber_params(n_nuc, diameter=diameter, rise=rise, nld=nld, ribs=ribs)
    # fiber_params[0,5] = 1.6
    fiber_dyad_frames = FMC.get_fiber_dyad_frames(fiber_params, nucl.d_index, nucl.nuc_type)
    FMC.cast_nucs_fiber(DNA, dyads, fiber_dyad_frames, cut=cut)


def MC(DNA, start=179, end=239, unwrap=-20, diameter=330.0, rise=100.0, nld=100, ribs=1, bert_factor = 1):

    #METROPOLIS-HASTINGS MONTE CARLO
    for i in range(start + unwrap, end - unwrap):

        Estart = FMC.Boltzmann_Energy(DNA.params[start + unwrap]) + FMC.Boltzmann_Energy(DNA.params[i])
        pars = DNA.params
        temp = (0.5 - np.random.rand()) * bert_factor * np.array([0.001, 0.001, 0.01, 0.001, 0.001, 0.01])
        pars[i] = pars[i]+temp
        # print temp
        newDNA = DNA.copy()
        newDNA.set_params(pars)
        GetCut(newDNA, nucl, 2, start - dyads[0] + unwrap, dyads, diameter=diameter, rise=rise, nld=nld, ribs=ribs)
        Enew = FMC.Boltzmann_Energy(newDNA.params[start + unwrap]) + FMC.Boltzmann_Energy(newDNA.params[i])

        if Enew - Estart < 0:
            DNA.set_params(newDNA.params)
            print Enew - Estart
        else:
            r = np.exp(-Enew/kBT)/np.exp(-Estart/kBT) #kBT is in Angstrom
            P = np.random.rand()

            if P>=r:
                DNA.set_params(newDNA.params)
                print Enew - Estart



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


def CurveDNA(apars, unwrap, DNA, nucl, start, end, n_nuc, dyads, diameter=330.0, rise=100.0, nld=100, ribs=1):
    """
    Curve the linker DNA
    """
    newpars = DNA.params

    newpars[start + unwrap:end - unwrap, :] = FMC.create_curved_linker(end - start - 2 * unwrap, apars)
    DNA.set_params(newpars)
    GetCut(DNA, nucl, n_nuc, start - dyads[0] + unwrap, dyads, diameter=diameter, rise=rise, nld=nld, ribs=ribs)


# ititiate program


# dna,dyads,nucl = FMC.create_fiber2(0, 2, 197,dna_file='simulations\\NLD\\20180220\\20180220_197_20_NLD.npz', nuc_file = '1KX5.3DNA', CL=False)
# FMC.plotdna2(dna,dyads,nucl)


min_unwrap = 0
max_unwrap = 10

date = time.strftime("%Y%m%d")
NRL = 197
file_name = 'simulations\\MC Metropolis-Hastings\\NLD series\\' + date + '_' + str(NRL)
if not os.path.exists(file_name):
    os.makedirs(file_name)


read_from_file = True

n_nuc = 2
n_bp = 147 + NRL
dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL)
# end = FMC.nuc_get_ends(dna, dyads[1], nucl)[0]
start = dyads[0] + (147 - nucl.d_index)
end = dyads[1] - nucl.d_index
# start = FMC.nuc_get_ends(dna, dyads[0], nucl)[1]

####Load or create a DNA sequence####
if read_from_file == True:

    ###Load DNA sequence from existing DNA pose###
    dna_file = 'simulations\\20180215_197_0_unwrap.npz'
    dna = HelixPose(params=np.load(dna_file)['params'], frame0=np.load(dna_file)['frame0'])

else:

    ###Create DNA sequence according to input parameters###
    # 197
    a = np.array([10.234, 0.312, 0.153, 0.973, 0.776])  # parameters in range we can probe with GUI
    a = np.array([ 11.6,     0.324,   3.625,   0.976,   0.795])  # on their way to be optimized...
#    a = np.array([11.617, 0.329, 3.672, 0.965, 0.800])  # 15/02/18 197 solenoid --> unwrap 0-10
    # a = np.loadtxt('simulations\\180215_197_0_unwrap.dat')

    unwrap = -5
    nld = 100

    link_len = end - start

    CurveDNA(a, unwrap, dna, nucl, start, end, n_nuc, dyads, nld = nld)
 
    
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
    print dna.params[start + unwrap]
    print Ecurr
 
    
    minimum = False

    # optimize the Bert/John parameters
    while minimum == False:
        a, Ecurr, minimum = FindMinimum(dna, nucl, a, dyads, start=start, end=end, unwrap=unwrap)

    Etemp = FMC.Boltzmann_Energy(dna.params[start + unwrap])
    print Etemp

    dna.write2disk(str(file_name) + '_' + str(nld) + '_nld_optimized')
    np.savetxt(str(file_name) + '_' + str(nld) + '_nld_input_params_optimized.dat', a)












##############Optimize linker DNA using Metropolis-Hastings Monte Carlo##############

nld_start = 100
nld_end = 17

n = 0
E_res = np.array([])
E_link_series = np.array([])

# for unwrap in range(min_unwrap,max_unwrap+1):
for nld in range(int(nld_start), int(nld_end-1),-1):
    n += 1
    
    # unwrap = -unwrap
    unwrap = 0

    # nld = 100

    E = []
    E_link = []
    E_bp = []

    j = 0

    print 'Start index of linker DNA: ' + str(start)
    print 'End index of linker DNA: ' + str(end)
    number_of_iterations = 1000
    while j < number_of_iterations:
        MC(dna, start, end, unwrap, nld=nld, bert_factor=1)

        # NLD
        print str(j+1) + ' of ' + str(number_of_iterations) + ' (' + str(nld) + ' nucleosome line density)'

        # unwrap
        # print str(j+1) + ' of ' + str(number_of_iterations) + ' (' + str(unwrap) + ' unwrapped base pairs)'

           # nice to include the energy calculated from the helixmc score-function
        Etemp = FMC.Boltzmann_Energy(dna.params[start + unwrap])
        print 'Energy of first base pair: ' + str(Etemp)
        E.append(Etemp)
        
        E_l = []
        for i in range(start + unwrap, end - unwrap):
            E_l.append(FMC.Boltzmann_Energy(dna.params[i]))
            
        E_link.append(np.sum(E_l))
        print 'Total energy of linker DNA: ' + str(np.sum(E_l))
        E_bp = np.concatenate((E_bp,np.array(E_l)))

        j += 1


# unwrap

#     np.savetxt(str(file_name) + '_' + str(unwrap) + '_unwrap_E.dat ', E)
#     np.savetxt(str(file_name) + '_' + str(unwrap) + '_unwrap_E_linkerDNA.dat ', E_link)
#     dna.write2disk(str(file_name) + '_' + str(unwrap) + '_unwrap')
#
#     E_bp = np.transpose(np.hsplit(E_bp, j))
#     np.savetxt(str(file_name) + '_basepair_energies_unwrap' + str(unwrap) + '.dat', E_bp)
#
#     E_res = np.concatenate((E_res, E))
#     E_link_series = np.concatenate((E_link_series, E_link))
#
# E_res = np.transpose(np.hsplit(E_res, n))
# np.savetxt(str(file_name) + '_energies_series.dat', E_res)
#
# E_link_series = np.transpose(np.hsplit(E_link_series, n))
# np.savetxt(str(file_name) + '_energies_linker_series.dat', E_link_series)



# NLD

    np.savetxt(str(file_name) + '_' + str(nld) + '_NLD_E.dat ', E)
    np.savetxt(str(file_name) + '_' + str(nld) + '_NLD_E_linkerDNA.dat ', E_link)
    dna.write2disk(str(file_name) + '_' +str(nld)+'_NLD')

    E_bp = np.transpose(np.hsplit(E_bp,j))
    np.savetxt(str(file_name) + '_basepair_energies_NLD' + str(nld) + '.dat',E_bp)

    E_res = np.concatenate((E_res, E))
    E_link_series = np.concatenate((E_link_series, E_link))

E_res = np.transpose(np.hsplit(E_res,n))
np.savetxt(str(file_name) + '_energies_series.dat',E_res)

E_link_series = np.transpose(np.hsplit(E_link_series,n))
np.savetxt(str(file_name) + '_energies_linker_series.dat',E_link_series)