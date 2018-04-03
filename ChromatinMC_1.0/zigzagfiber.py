# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:08:55 2016

@author: Visscher
"""

import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import NucleosomeMC as nuc
import FiberMC as FMC
from helixmc.pose import HelixPose
from helixmc import util


n_nuc = 8
NRL = 167
n_bp = NRL * n_nuc
dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL, unwrap = None, compute_tw_wr = False, nuc_file = '1KX5.3DNA')
tf = nucl.com_d_transform

def ribbon(r, cg, theta0, N, ld, h0):
    '''
    Creates dyad frames for nucleosomes in a ribbon    
    
    Parameters 
    ----------
    r : float 
        fiber radius
    theta0 : float
        starting angle
    N : int
        number of nucleosomes
    ld : float
        distance between two nucleosomes in the same ribbon
    h0: float
        beginning height of the spiral
    
    Returns
    -------
    coords : list
        list of nucleosome coords
    '''
    
    
    def theta(s):
        return s/r + theta0
    
    def znuc(s):
        return 1.0/np.sqrt(1+cg**2)*np.array([-np.sin(theta(s)),\
                np.cos(theta(s)), cg])
                
    def xnuc(s): 
        return np.array([-np.cos(theta(s)), -np.sin(theta(s)), 0])
        
    def ynuc(s):
        return np.cross(znuc(s), xnuc(s))
        
    def fnuc(s):
        return np.array([xnuc(s),ynuc(s),znuc(s)])
    
    def onuc(s):
        x = r*np.cos(theta(s))
        y = r*np.sin(theta(s))
        z = s*cg + h0
        return np.array([x,y,z])
        
    ss = np.arange(0,N*ld, ld) + h0/cg
    onucs = []
    fnucs = []
    cnucs = []
    for s in ss:
        o = onuc(s)
        f = fnuc(s)
        onucs.append(o)
        fnucs.append(f)
        cnucs.append(nuc.of2coords(o,f))
    

    dfs = []
    for i in range(len(onucs)):


       
        dna_params, dna_o, dna_f0, c = FMC.GetNuc(norigin = onucs[i], nframe = fnucs[i])
        n_dna = HelixPose(dna_params, frame0 = dna_f0)

        dyad_coords = nuc.of2coords(n_dna.coord[nucl.d_index]+dna_o, np.transpose(n_dna.frames[nucl.d_index]))
        dfs.append(dyad_coords)
        

    return dfs, cnucs

def connect_nucs(N, Nrib, Nstep):
    '''
    Creates a sequence of nucleosome ribbon number and height indices
    
    Parameters
    ----------
    N : int
        Total number of nucleosomes in the fiber
    Nrib : int
        Number of ribbons in the fiber
    Nstep : int
        Number of steps 
    
    Returns
    -------
    nuc_list : list
        sequence of nucleosome ribbon number and height indices
    '''
    i = 0 #ribbon number
    j = 0 #index in ribbon
    
    nuc_list = []
    
    totalnucs = 0
    while totalnucs < N:
        
        if totalnucs % Nrib == 0 and totalnucs != 0:
            j += 1
        nuc_list.append([i,j])
        i = (i + Nstep) % Nrib
        totalnucs += 1
    
    return nuc_list
    
    
    

def fiber(D, Nrib, Nstep, Nnuc):
    '''
    Takes in fiber parameters and constructs a list of dyad frames based on 
    these parameters
    
    Parameters
    ----------
    D : float
        diameter of the fiber
    Nrib : int
        number of ribbons
    Nstep : int
        number of nucleosomes between two linked nucleosomes
    Nnuc : int
        total number of nucleosomes, multiple of Nrib
    rise : float
        rise per nucleosome
    NLD : float
        nucleosome line distance
    
    Returns
    -------
    dyad_frames : list
        a list of dyad_frames
    '''
    a = 115.0 #diameter of NCP in angstrom
    b = 60.0 #height of NCP in angstrom
    r = D/2.0
    
    
    theta0s = np.arange(0, 2*np.pi, 2*np.pi/Nrib)
    alpha = 2*b/(D-a)*(1-(a*Nrib/(np.pi*(D-a)))**2)
    hs = np.arange(0,Nrib*b*np.sin(alpha),b*np.sin(alpha))
    cg = a*Nrib/(np.pi*(D-a))
    spiralcs = []
    for i in range(Nrib):
        spiralcs.append(ribbon(r, cg, theta0s[i], Nnuc/Nrib, b, hs[i]))
        
        
    
    nuc_list = connect_nucs(Nnuc, Nrib, Nstep)
    dyad_frames = []
    for lst in nuc_list:
        i = lst[0]
        j = lst[1]
        dyad_frames.append(spiralcs[i][j])
    

    

        
    
    
    return dyad_frames
    
def fiber2(D, Nnuc, rise, NLD):
    '''
    Takes in fiber parameters and constructs a list of dyad frames based on 
    these parameters
    
    Parameters
    ----------
    D : float
        diameter of the fiber
    Nrib : int
        number of ribbons
    Nstep : int
        number of nucleosomes between two linked nucleosomes
    Nnuc : int
        total number of nucleosomes, multiple of Nrib
    rise : float
        rise per nucleosome
    NLD : float
        nucleosome line distance
    
    Returns
    -------
    dyad_frames : list
        a list of dyad_frames
    '''
    Nrib = 2
    Nstep = 1

    r = D/2.0
    theta0s = np.array([0,np.pi])
    hs = np.array([0, NLD])
    cg = 2.0*NLD/rise
    
    
    spiralcs = []
    nucs = []
    for i in range(Nrib):
        spiralcs.append(ribbon(r, cg, theta0s[i], Nnuc/Nrib, rise, hs[i])[0])
        nucs.append(ribbon(r, cg, theta0s[i], Nnuc/Nrib, rise, hs[i])[1])
        
    
    nuc_list = connect_nucs(Nnuc, Nrib, Nstep)
    dyad_frames = []
    nuc_frames = []
    for lst in nuc_list:
        i = lst[0]
        j = lst[1]

    
        dyad_frames.append(spiralcs[i][j])
        nuc_frames.append(nucs[i][j])

    

        
    
    
    return dyad_frames, nuc_frames
    


def plotdna(DNA, dfs):
    '''
    Save DNA frames using matplotlib + mplot3d

    Parameters
    ----------
    plt_range : int, optional
        xyz range (nm), default is DNA length
    force : float, optional 
        the force applied to the strand of DNA, if not None will add it to the 
        title of the figure
    n_MC : number of the MC cycle. If not None will be added to the title figure

      '''
    plt_range = 40
    plt.close()
    rb_width=2.0
    rb_width=2.0
    color='kb' 
    
    rb_vec = DNA.rb_vec
    coord = DNA.coord/10
    frames = np.transpose(DNA.frames, axes = (0,2,1))
    tf2 = nuc.get_transformation(dfs[0], nuc.of2coords(DNA.coord[83], np.transpose(DNA.frames[83])))
    for i in range(len(dfs)):
        dfs[i] = nuc.apply_transf_coords(dfs[i],tf2)
    
    bb1 = coord - rb_vec * rb_width * 0.5
    bb2 = coord + rb_vec * rb_width * 0.5
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.plot3D(bb1[:, 0], bb1[:, 1], bb1[:, 2], color[0] + '-')
    ax.plot3D(bb2[:, 0], bb2[:, 1], bb2[:, 2], color[0] + '-')

    

    
            




#    for i in xrange(bb1.shape[0]):
#        ax.plot3D(
#            np.array([bb1[i, 0], bb2[i, 0]]),
#            np.array([bb1[i, 1], bb2[i, 1]]),
#            np.array([bb1[i, 2], bb2[i, 2]]),
#            color[1] + '-')
            
        
    color = ['b','g','r']
    for j in range(len(dfs)):
        dfs[j] = dfs[j]/10
        for i in np.arange(1,4):
            ax.plot3D(
                [dfs[j][0,0],dfs[j][0,0]+(dfs[j][i,0]-dfs[j][0,0])*30],
                [dfs[j][0,1],dfs[j][0,1]+(dfs[j][i,1]-dfs[j][0,1])*30],
                [dfs[j][0,2],dfs[j][0,2]+(dfs[j][i,2]-dfs[j][0,2])*30], color[i-1])
                
    for d in dyads:
        for i in np.arange(0,3):
            ax.plot3D(
                [coord[d,0], coord[d,0] + 3*frames[d,i,0]],
                [coord[d,1], coord[d,1] + 3*frames[d,i,1]],
                [coord[d,2], coord[d,2] + 3*frames[d,i,2]], color[i-1])
    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')

    if plt_range is None:
        plt_range = DNA.n_bp*0.34      
    ax.set_xlim3d(-plt_range/2,plt_range/2)
    ax.set_ylim3d(-plt_range/2,plt_range/2)
    ax.set_zlim3d(0,plt_range)
    
    
    
def plotnucl(nucl):
    '''
    Save DNA frames using matplotlib + mplot3d

    Parameters
    ----------
    plt_range : int, optional
        xyz range (nm), default is DNA length
    force : float, optional 
        the force applied to the strand of DNA, if not None will add it to the 
        title of the figure
    n_MC : number of the MC cycle. If not None will be added to the title figure

      '''
    plt_range = 35
    plt.close()
    rb_width=2.0
    rb_width=2.0
    color='kb' 
    
    DNA = nucl.DNApose
    rb_vec = DNA.rb_vec
    coord = DNA.coord/10
    
    color = ['b','g','r']    
    
    bb1 = coord - rb_vec * rb_width * 0.5
    bb2 = coord + rb_vec * rb_width * 0.5
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.plot3D(bb1[:, 0], bb1[:, 1], bb1[:, 2], color[0] + '-')
    ax.plot3D(bb2[:, 0], bb2[:, 1], bb2[:, 2], color[0] + '-')

    n_coords = nucl.n_coords/10.0

    
    for i in np.arange(1,4):
        ax.plot3D(
            [n_coords[0,0],n_coords[0,0]+(n_coords[i,0]-n_coords[0,0])*30],
            [n_coords[0,1],n_coords[0,1]+(n_coords[i,1]-n_coords[0,1])*30],
            [n_coords[0,2],n_coords[0,2]+(n_coords[i,2]-n_coords[0,2])*30], color[i-1])
            
    d_coords = nucl.d_coords/10.0

    
    for i in np.arange(1,4):
        ax.plot3D(
            [d_coords[0,0],d_coords[0,0]+(d_coords[i,0]-d_coords[0,0])*30],
            [d_coords[0,1],d_coords[0,1]+(d_coords[i,1]-d_coords[0,1])*30],
            [d_coords[0,2],d_coords[0,2]+(d_coords[i,2]-d_coords[0,2])*30], color[i-1])

    for i in xrange(bb1.shape[0]):
        ax.plot3D(
            np.array([bb1[i, 0], bb2[i, 0]]),
            np.array([bb1[i, 1], bb2[i, 1]]),
            np.array([bb1[i, 2], bb2[i, 2]]),
            color[1] + '-')

    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')

    if plt_range is None:
        plt_range = DNA.n_bp*0.34      
    ax.set_xlim3d(-plt_range/2,plt_range/2)
    ax.set_ylim3d(-plt_range/2,plt_range/2)
    ax.set_zlim3d(0,plt_range)


dfs, nfs = fiber2(250.0, n_nuc, 100, 15)
nfs = np.array(nfs)
##print dfs

#plt_range = 350
#plt.close()
#fig = plt.figure()
#ax = p3.Axes3D(fig)
#
#color = 'rbg'
#for j in range(len(nfs)):
#    dfs[j] = dfs[j]/10
#    for i in np.arange(1,4):
#        ax.plot3D(
#            [nfs[j][0,0],nfs[j][0,0]+(nfs[j][i,0]-nfs[j][0,0])*30],
#            [nfs[j][0,1],nfs[j][0,1]+(nfs[j][i,1]-nfs[j][0,1])*30],
#            [nfs[j][0,2],nfs[j][0,2]+(nfs[j][i,2]-nfs[j][0,2])*30], color[i-1])
#
#ax.set_xlim3d(-plt_range/2,plt_range/2)
#ax.set_ylim3d(-plt_range/2,plt_range/2)
#ax.set_zlim3d(0,plt_range)
##print dyads
#FMC.cast_nucs_fiber(dna, dyads, dfs, cut = None)
#
#
#plotdna(dna, nfs)

