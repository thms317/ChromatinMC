# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:55:46 2016

@author: Visscher
"""

import numpy as np

import FiberMC as FMC

from Cscore import ScoreTweezers
from multiprocessing import Pool, freeze_support
import time
import os
from helixmc import random, kBT

#define fiber
n_nuc = 1
NRL = 197
n_bp = NRL*n_nuc + 500
#define simulation
steps = 500

def write2file(output, header = ['force (pN)', 'z (nm)']):
    '''
    Write the results of a simulated experiment to a text file
    
    Parameters
    ----------
    output : array 
        the output to write to the file
    header : list of str
        the header of the file
        
    '''
    datfile = FMC.get_filename(ext = 'dat', incr = True)
    FMC.write_dat(output, file = datfile, header = header)

def SingleMC(force):
    '''
    Run a MC simulation for a single datapoint
    
    Parameters
    ----------
    force : float
        the force 
    '''
    dna, dyads, nucl = FMC.create_fiber(0, 4, 197, 'DinucleosomePoses\\1KX5NRL197.npz')
    n_bp = len(dna.coord)
    print n_bp
    unwrap = -30
    scorefxn = ScoreTweezers(force)
    best_is = {}
    for d in dyads:
        best_is[d] = 6
    bps = np.arange(dyads[0]-(75+ unwrap))
    g = 10*kBT
    fps = []
    print FMC.MakeFiberHelix(dna, dyads, best_is)
    for i in range(len(dyads)-1):
        bps = np.concatenate((bps, np.arange(dyads[i]+ (71+unwrap), dyads[i+1] -(75+ unwrap))))
    bps = np.concatenate((bps, np.arange(dyads[-1]+(71+unwrap), n_bp -1)))
    print bps
    bps = np.arange(n_bp -1)
    for i in bps:
        
        FMC.MC_move(i, scorefxn, dna, dyads, nucl,best_is)
    LK0 = dna.link_exact
    dLK = []
    z = []
    uw = []
    nucl.fixed_g = g

    for j in range(steps):
        acc = 0
        runtime = time.time()
        for i in bps:
            
            
            acc += FMC.MC_move(i, scorefxn, dna, dyads, nucl, best_is)
        z.append(dna.coord_terminal[-1]/10)
        uw.append(FMC.count_unwrap2(dna, dyads, nucl, best_is[dyads[0]]))
        dLK.append(dna.link_exact-LK0)
        print 'Step %d force %.2f ' %(j,force)
        print time.time() - runtime
        #print FMC.MakeFiberHelix(dna, dyads, best_is)
        print 'accepted %.2f procent of moves' %(100.0*float(acc)/len(bps))
        dna.write2disk('%.0ftestrun%d'%(force*100,j))
        fps.append(FMC.MakeFiberHelix(dna, dyads, best_is))
        if j%25 == 0:
            fps = np.array(fps)
            np.save('fps%d' %(j/25), fps)
            fps = []
            z = np.array(z)
            np.save('z%d' %(j/25), z)
            z = []
    fps = np.array(fps)
    z = np.array(z)
    uw = np.array(uw)
    dLK = np.array(dLK)
#    np.save('force%d'%(force*100),z)
#    np.save('uw%dKT%d'%(g/kBT,force) ,uw)
#    np.save('dLK%dKT%d'%(g/kBT,force) ,dLK)
    np.save('fpstest', fps)
    output = np.array([force,np.mean(z),np.std(z)])
    #np.savetxt('out%.0f.npy' %force*100, np.array(output))
    write2file(output)
    #dna.write2disk('NRL1701KX5render')
    return output
    
def MultiMC(k):
    '''
    Run a MC simulation for multiple datapoints
    
    Parameters
    ----------
    forces : iterable of floats
        the forces at which to do MCs
        
    '''
    force1 = np.arange(0,5)
    force2 = np.arange(5,10)
    force3 = np.arange(10,15)
    force4 = np.arange(15,20)


    
    if k == 1:
        forces = force1
        dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL, unwrap = 0, nuc_file = '1KX5.3DNA')
        nucl.fixed_g = 10*kBT
        g = 3
    if k == 2:
        forces = force2
        dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL, unwrap = 0, nuc_file = '1KX5.3DNA')
        nucl.fixed_g = 10*kBT
        g = 10
    if k == 3:
        forces = force3
        dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL, unwrap = 0, nuc_file = '1KX5.3DNA')
        nucl.fixed_g = 10*kBT

    if k == 4:
        forces = force4
        dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL, unwrap = 0, nuc_file = '1KX5.3DNA')
        nucl.fixed_g = 10*kBT
        

    best_is = {}
    datfile = FMC.get_filename(file = 'MultiMC' + str(k), ext = 'dat', incr = True)
    for dyad in dyads:
        best_is[dyad] = 6
    
    for f in forces:
        scorefxn = ScoreTweezers(f)
        uw = []
        dLK = []
        for i in range(n_bp-1):
            #fiberparams = FMC.MakeFiberHelix(dna, dyads)
            FMC.MC_move(i, scorefxn, dna, dyads, nucl, best_is)
    
        z = []
        
    
        for j in range(steps):
            runtime = time.time()
            for i in range(n_bp-1):
                
                #fiberparams = FMC.MakeFiberHelix(dna, dyads)
                FMC.MC_move(i, scorefxn, dna, dyads, nucl,best_is)
            z.append(dna.coord_terminal[-1]/10)
            print 'Step %d force %.2f ' %(j,f)
            print time.time() - runtime
            uw.append(FMC.count_unwrap2(dna, dyads, nucl, best_is[dyads[0]]))
            z.append(dna.coord_terminal[-1]/10)
            
        z = np.array(z)
        uw = np.array(uw)
        dLK = np.array(dLK)

#        np.save('uw%dKT%d'%(g,f) ,uw)
#        np.save('dLK%dKT%d'%(g,f) ,dLK)

        result = np.array([f, np.mean(z), np.std(z)])
        FMC.write_dat(result, file = datfile, header = ['force (pN)', 'mean[z] (nm)', 'std[z] (nm)'])
    
    #suw = np.array(uw)
   
    return
    
    
def NucBreathing(force):
    '''
    perform a nucleosome breathing simulation
    '''
    scorefxn = ScoreTweezers(force)
    dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL, unwrap = 0)
    

#force4 = 7 - np.linspace(0,4,12)
#force3 = 3 + np.linspace(0,4,12)
#force1 = np.linspace(0,1,12)
#force2 = np.linspace(1,2,12)
#
#forces = [force1, force2, force3, force4]
#for i in range(0,5):
#    SingleMC(float(i), 3*kBT)
#    SingleMC(float(i), 10*kBT)
#SingleMC(float(5), 3*kBT)

if __name__ == '__main__':
    freeze_support()
    pool = Pool(4)
    i = [1,2,3,4]
    
    

    
    runtime1 = time.time()
    print runtime1
    result1 = pool.map(MultiMC, i)
    runtime1 = time.time() - runtime1
    print runtime1
  


