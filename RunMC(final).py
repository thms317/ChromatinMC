# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 10:37:51 2016

@author: Visscher
"""
from multiprocessing import Pool, freeze_support
import os
import numpy as np
from Cscore import ScoreTweezers
import time
from datetime import datetime
from helixmc import random, kBT

import FiberMC as FMC
import NucleosomeMC as NMC

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
    
    
def write_in_folder(filename, data):
    '''
    Write a numpy array in a data folder and create a data folder if it does not
    exist yet
    
    Parameters
    ----------
    filename : string
        the name of the file to write
    data : ndarray
        the array to be saved
        
    Returns
    -------
    
    '''
    #check if folder exists
    i = datetime.now()
    currdir = os.getcwd()
    
    folder = '\\data' + i.strftime('%m%d') + '\\'
    d = currdir + folder
    print(d)
    if not os.path.exists(d):
        os.makedirs(d)
    np.save(d+filename, data)

def create_bp_seq(dyads, best_is, n_bp, nucl):
    '''
    Create the sequence in which to run the Monte Carlo algorithm
    
    Parameters
    ----------
    dyads : list
        list of dyad basepair indices
    best_is : dictionary
        dictionary of the mean wrapped fixed base pair per nucleosome
        dyad indices function as keys
    n_bp : int
        number of basepairs        
        
    Returns
    -------
    bps : ndarray
        the sequence of indices to run the MC algorithm
    
    '''
    
    bps = np.array([])
    start = 0
    for d in dyads:
        end = d + nucl.fixed_i[-1] + 2
        bps = np.concatenate((bps, np.arange(start,nucl.fixed_i[best_is[d]]+d), \
            np.arange(nucl.fixed_i[best_is[d]]+d,d + nucl.fixed_i[-1] + 2)[::-1]))
        start = end
    
    bps = np.concatenate((bps, np.arange(start, n_bp)))
    
    return bps

 
    
    

def Run_Force_Extension(forces, steps, n_nuc = 1, n_bp = None, NRL = 197, g = 10, nuc_type = '1KX5', unwrap = 0):
    '''
    Run a force extension simulation
    
    Parameters
    ----------
    forces : list
        list of the forces at which to simulate
    n_nuc : int
        number of nucleosomes
    n_bp : int
        number of base pairs
    g : float
        Value of the adsorbion energy in kBT
    nuc_type : string
        The type of nucleosome crystal structure
    unwrap : int
        Number of base pairs to be unwrapped at the start
        (in case of mononucleosome pulling)
    '''
    if n_bp == None:
        n_bp = (n_nuc - 1)*NRL + 167
    
    if n_nuc == 1:
        dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL, unwrap = unwrap, nuc_file = nuc_type + '.3DNA')

    else:
        dna, dyads, nucl = FMC.create_fiber(n_bp -  (n_nuc - 1)*NRL + 147, \
        n_nuc, NRL, 'DinucleosomePoses\\'+ nuc_type + 'NRL%d.npz' %NRL)

    nucl.fixed_g = g*kBT
    
    best_is = {}
    for d in dyads:
        best_is[d] = 6

    datfile = FMC.get_filename(file = 'Force_extension', ext = 'dat', incr = True)   
    for f in forces:
        scorefxn = ScoreTweezers(f)      
        Lk = []
        Wr = []
        z = []    
        for j in range(steps):
            start_time = time.time()
            bps = create_bp_seq(dyads, best_is, n_bp, nucl)
            acc = 0
            for i in bps:
                acc += FMC.MC_move(i, scorefxn, dna, dyads, nucl, best_is)
            
            Lk.append(dna.link_exact/(2*np.pi))
            Wr.append(dna.writhe_exact/(2*np.pi))
            z.append(dna.coord_terminal[-1]/10)
            print('Step %d of %d at force %.2f completed!' %(j+1, steps, f))
            print('Required time : %.2f seconds' %(time.time()-start_time))
            print('Accepted %.2f percent of moves' %(float(acc)/n_bp*100))
    
        Lk = np.mean(np.array(Lk))
        Lkerr = np.std(np.array(Lk))        
        Wr = np.mean(np.array(Wr))
        Wrerr = np.std(np.array(Wr))
        z = np.mean(np.array(z))
        zerr = np.std(np.array(z))
        result = np.array([f, z, zerr, Lk, Lkerr, Wr, Wrerr])
        FMC.write_dat(result, file = datfile, header = ['force (pN)', 'mean[z] (nm)', 'std[z] (nm)', 'lk', 'lkerr', 'wr', 'wrerr'])


def Nuc_Breathing(force, steps, g, n_bp = 647, nuc_type = '1KX5', unwrap = 0):
    '''
    Parameters
    ----------
    force : int
        The force at which to run the simulation

    n_bp : int
        number of base pairs
    g : float
        Value of the adsorbion energy in kBT
    nuc_type : string
        The type of nucleosome crystal structure
    unwrap : int
        Number of base pairs to be unwrapped at the start
        (in case of mononucleosome pulling)
    '''
    dna, dyads, nucl = FMC.create_dna(n_bp, 1, 197, nuc_file = nuc_type + '.3DNA', unwrap = unwrap)

    best_is = {dyads[0]:6}
    scorefxn = ScoreTweezers(force)
    nucl.fixed_g = g*kBT
    
    uw = []    
    
    for j in range(steps):
        start_time = time.time()
        bps = create_bp_seq(dyads, best_is, n_bp, nucl)
        acc = 0
        for i in bps:
            acc +=FMC.MC_move(i, scorefxn, dna, dyads, nucl, best_is)
        uw.append(FMC.count_unwrap2(dna, dyads, nucl, best_is[dyads[0]]))
        print('Step %d of %d at force %.2f completed!' %(j+1, steps, force))
        print('Required time : %.2f seconds' %(time.time()-start_time))
        print('Accepted %.2f percent of moves' %(float(acc)/steps*100))
    uw = np.array(uw)
    return uw      
         
        
        
if __name__ == '__main__':
        Run_Force_Extension(np.arange(0,11,1), 10, n_nuc = 1, n_bp = None, NRL = 197, g = 10, nuc_type = '1KX5', unwrap = 0)
        
        
        
        
        
        
        
        
        
        
        
    

