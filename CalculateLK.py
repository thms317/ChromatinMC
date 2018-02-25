# -*- coding: utf-8 -*-
"""
Created on Wed May 04 11:46:27 2016

@author: Visscher
"""
from helixmc.pose import HelixPose
import numpy as np
import NucleosomeMC as NMC
import FiberMC as FMC
from matplotlib import pyplot as plt

def GetPose(filename, drc = ''):
    '''
    Load a dinucleome HelixPose object from .npz file
    
    Parameters
    ----------
    filename : string
        name of the .npz file
    drc : string
        directory in which the file can be found (if it is not in the same 
        directory as this script)
    
    Returns
    -------
    DNA : HelixPose object
        the HelixPose of a dinucleosome
    dyads : ndarray
        numpy array with dyad indices
    nucl : NucPose object
        the NucPose of the used nucleosome
    '''
    fname = drc + '\\' + filename
    
    NRL = int(filename.split('NRL')[1].split('h')[0])
    nuctype = filename.split('NRL')[0]
    #h = int(filename.split('NRL')[1].split('h')[1].split('.')[0])
    
    DNA = HelixPose(params = np.load(fname)['params'] , frame0 = np.load(fname)['frame0'])
    n_bp = len(DNA.params)
    dyads = np.asarray(NRL*(np.arange(0, 2, 1)-(2-1)/2.0))
    dyads = (dyads+n_bp/2).astype(int)
    nucl = NMC.NucPose()
    nucl.from_file(nuctype+'.3DNA')
    return DNA, dyads, nucl
    
   
def SplitParams(DNA, dyads, nucl):
    '''
    Split the dinucleosome pose in a linker and a nucleosome part
    
    Parameters
    ----------
    DNA : HelixPose object
    
    dyads : ndarray
    
    nucl : NucPose object
    
    Returns
    -------
    n_params : ndarray
        array of params of the nucleosomal dna
    l_params : ndarray
        array of params of the linker DNA
    '''
    n_params = DNA.params[dyads[0]-nucl.d_index:dyads[0]-nucl.d_index+145]    
    l_params = DNA.params[dyads[0]-nucl.d_index+146: dyads[1]+nucl.d_index]

    return n_params, l_params


#DNA, dyads, nucl = GetPose('1KX5NRL197h25.npz', 'DinucleosomePoses')
#
#
#
#
#
#
#
#n_pars, l_pars = SplitParams(DNA, dyads, nucl)
#DNA_nuc = DNA.copy()
#DNA_link = DNA.copy()
#DNA_test = DNA.copy()
#DNA_nuc.set_params(n_pars)
#DNA_link.set_params(l_pars)
#DNA_test.set_params(np.tile(np.concatenate((n_pars,l_pars)), (6,1)))
#
##print len(n_pars)
##print len(l_pars)
#
#
#
#print DNA_test.link_exact/(2*np.pi*12) - 197.0/10.4
#print 197.0/10.4


NRL = 167
n_nuc = np.arange(2,100)
lk = []
wr = []
for n in n_nuc:
    dna, dyads, nucl = FMC.create_fiber(0, n, NRL, dna_file = 'DinucleosomePoses\\4qlcNRL%d.npz'%NRL)
    newpars = dna.params[dyads[0]:dyads[-1]]    
    print n
    lk.append(dna.link_exact/(2*n*np.pi))
    wr.append(dna.writhe_exact/(2*n*np.pi))

lk = np.array(lk)
wr = np.array(wr)
np.save('4qlc167lk2',lk)
np.save('4qlc167wr2',wr)

NRL = 197
n_nuc = np.arange(2,100)
lk = []
wr = []
for n in n_nuc:
    dna, dyads, nucl = FMC.create_fiber(0, n, NRL, dna_file = 'DinucleosomePoses\\4qlcNRL%d.npz'%NRL)
    newpars = dna.params[dyads[0]:dyads[-1]]    
    print n
    lk.append(dna.link_exact/(2*n*np.pi))
    wr.append(dna.writhe_exact/(2*n*np.pi))

lk = np.array(lk)
wr = np.array(wr)
np.save('4qlc197lk2',lk)
np.save('4qlc197wr2',wr)

NRL = 170
n_nuc = np.arange(2,100)
lk = []
wr = []
for n in n_nuc:
    dna, dyads, nucl = FMC.create_fiber(0, n, NRL, dna_file = 'DinucleosomePoses\\1KX5NRL170.npz')
    newpars = dna.params[dyads[0]:dyads[-1]]    
    print n
    lk.append(dna.link_exact/(2*n*np.pi))
    wr.append(dna.writhe_exact/(2*n*np.pi))

lk = np.array(lk)
wr = np.array(wr)
np.save('1KX5170lk2',lk)
np.save('1KX5170wr2',wr)

##for n in np.arange(23):
##    print n+2, lsrad[n], lspnrad[n]
#np.save(FMC.nuc_name.split('.')[0]+'NRL%dlink'%NRL,lspnrad)  
#
#NRL = 167
#n_nuc = np.arange(2,46)
#for n in n_nuc:
#    dna, dyads, nucl = FMC.create_fiber(0, n, NRL, dna_file = 'DinucleosomePoses\\'+ FMC.nuc_name.split('.')[0]+'NRL%d.npz'%NRL)
#    print n
#    ls[n-2] = dna.link_exact
#    lspn[n-2] = ls[n-2]/n
#    lsrad[n-2] = ls[n-2]/(np.pi*2)
#    lspnrad[n-2] = lspn[n-2]/(np.pi*2)

#for n in np.arange(23):
#    print n+2, lsrad[n], lspnrad[n]
#np.save(FMC.nuc_name.split('.')[0]+'NRL%dlink'%NRL,lspnrad)      
#plt.plot(n_nuc,lspnrad)
#plt.xlabel('Fiber length (# of nucleosomes)')
#plt.ylabel('Linking Number')
#plt.title('Linking Number')
#plt.show()    
#i = 100
#
#dyads =  FMC.CreateFiber(DNA, dyads[0], dyads[1], n = i)
#n_bp = len(DNA.params) - 25
#exp_twist = float(n_bp)/10.4
#tot_twist = 0.5*DNA.link_exact/np.pi - 25/10.4
#
#print (tot_twist - exp_twist)/float(i) + 1.7/2.0

#names = ['1KX5NRL197.npy', '1KX5NRL170.npy', '4qlcNRL197.npy', '4qlcNRL167.npy', \
#'4qlcNRL197link.npy', '4qlcNRL167link.npy','1KX5NRL197link.npy', \
#'1KX5NRL170link.npy']
#
#
#for n in names[0:4]:
#    plt.plot(n_nuc,np.load(n), 'bo')
#    plt.xlabel('Fiber length (# of nucleosomes)')
#    plt.ylabel('Writhe')
#    plt.title('Writhe (' + n.split('.')[0] + ')')
#    plt.savefig(n.split('.')[0])
#    plt.close()
#
#for n in names[4:]:
#    plt.plot(n_nuc,np.load(n), 'bo')
#    plt.xlabel('Fiber length (# of nucleosomes)')
#    plt.ylabel('Linking Number')
#    plt.title('Linking Number (' + n.split('.')[0]+')')
#    plt.savefig(n.split('.')[0])
#    plt.close()



