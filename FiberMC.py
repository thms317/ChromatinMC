# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 14:03:53 2015

@author: noort

edit: bert

Add functionalities for chromatin fibers to HelixMC


"""
import numpy as np
from helixmc.pose import HelixPose
import NucleosomeMC as nuc
import mpl_toolkits.mplot3d.axes3d as p3
from helixmc import util
from Cscore import ScoreTweezers
from helixmc.random_step import RandomStepSimple
import matplotlib.pyplot as plt
import os, sys
from helixmc import random, kBT
import math


nuc_name = '1KX5.3DNA'
nuc_dum = nuc.NucPose()
nuc_dum.from_file(nuc_name)
nuc_d = nuc_dum.d_index
random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')
fixed = nuc_dum.fixed_i

###
###             UTILITIES AND FRAME/DATA SAVING
###

def make_gif(files, out_file, delay=20, repeat=True):
    """
    Uses imageMagick to produce an animated .gif from a list of
    picture files.
    """
    loop = -1 if repeat else 0
    imagefiles = files[0].split('_')[0]+'*.'+files[0].split('.')[-1] 
    command = 'convert -delay %d -loop %d %s %s' % (delay, loop, imagefiles, out_file)
    os.system(command)
    if os.path.isfile(out_file):
        for f in files:
            os.remove(f)
    else:
        print 'Error: ImageMagick did not create gif using command:'
        print command
        
def get_directory():    
    """
    Get or create default directory
    """
    from datetime import datetime

    username = os.environ.get( "USERNAME" )
    i = datetime.now()
    d = 'C:\\Users\\'+ username + '\\data\\' + i.strftime('%Y%m%d') + '\\' 
    d = os.path.dirname(d)
    if not os.path.exists(d):
        os.makedirs(d)
    return d
        
def get_filename(file = None, ext = 'dat', incr = False):    
    """
    Create default filename in default directory
    """
    from os.path import isfile, join

    
    d = get_directory()
    if file is None:
        base='data'
    else:
        base=file.split('.')[0]
        base=base.split('_')[0]

    allfiles = [f for f in os.listdir(d) if isfile(join(d, f))]
    cur_num = 0
    for f in allfiles:
        if f.split('_')[0] == base:
            dot = f.find('.')
            start = f.find('_')+1
            cur_num = np.max([cur_num, int(f[start:dot]) ])
    if incr: cur_num += 1

    out = d+"\\"+base+"_%04d." % cur_num  + ext 
    return out
    
def write_dat(dat, file = None, header = None):
    """
    Write tab seperated txt file, including headers for each column
    """
    if file is None:
        file = get_filename(ext = 'dat', incr = True)
    if not os.path.isfile(file):
        with open(file, 'a') as the_file:
            the_file.write(('\t'.join([str(x) for x in header]))+'\n')
    with open(file,'a') as f_handle:
        np.savetxt(f_handle, dat[None], delimiter="\t", fmt = '%.3f' )
    return
    
def SaveDNAFrame(DNA, files=None, plt_range=None, coord2 = [], coord3 = [], force = None, n_MC = None):
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
    plt.close()
    rb_width=2.0
    rb_width=2.0
    color='kb' 
    rb_vec = DNA.rb_vec
    coord = DNA.coord/10
      
    bb1 = coord - rb_vec * rb_width * 0.5
    bb2 = coord + rb_vec * rb_width * 0.5
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.plot3D(bb1[:, 0], bb1[:, 1], bb1[:, 2], color[0] + '-')
    ax.plot3D(bb2[:, 0], bb2[:, 1], bb2[:, 2], color[0] + '-')

    if len(coord2) != 0:
        coord2 = coord2/10
    for i in coord2:
        ax.plot3D(i[:,0], i[:,1], i[:,2], 'ro', ms=10)
    if len(coord3) != 0:
        coord3 = coord3/10
    for i in coord3:
        ax.plot3D(i[:,0], i[:,1], i[:,2], 'go', ms =10)
#        print i[:,0], i[:,1], i[:,2]

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
    figtitle = ""
    maketitle = False
    if force is not None:
        figtitle += "Force : "+"{:.2f}".format(force) + " pN "
        maketitle = True
    if n_MC is not None:
        figtitle += "Number of MC-cycle : " +"{:.0f}".format(n_MC)
        maketitle = True
    if maketitle:
        plt.title(figtitle)
        
    new_file = get_filename(file='tmp', ext = 'jpeg', incr = 1)
    ax.figure.savefig(new_file)
    if files is None:
        files = [new_file]
    else:
        files.append(new_file)

    return files
    

###
###             DNA AND NUCLEOSOME CREATION
###

def create_dna(n_bp, n_nucs, NRL, unwrap = None, compute_tw_wr = False, nuc_file = nuc_name):
    '''
    Create random DNA configuration, including nucleosomes

    Parameters
    ----------
    n_bp : number of basepairs
    n_nucs : number of nucleosomes
    NRL : nucleosome repeat length
    unwrap : number of basepaires wrapped on the nucleosome. Negative numbers
        define number of unwrapped basepairs on each side
    
    Returns
    ----------
    DNA : instance of HelixMC
    dyads : bp position of dyads
    nucl : instance of NucleosomeMC
    '''
    dyads = np.asarray(NRL*(np.arange(0, n_nucs, 1)-(n_nucs-1)/2.0))
    dyads = (dyads+n_bp/2).astype(int)
    
    random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')
    params = np.tile(random_step.params_avg, (n_bp, 1))
    
    # Get nucleosome
    nucl=nuc.NucPose()
    nucl.from_file(nuc_file)

    
    # Insert nucleosomes in DNA
    params, free_dna = insert_nucs(params, nucl.params, dna_dyads=dyads, nuc_dyad = nucl.d_index, wrapped_bp = unwrap)
    DNA = HelixPose(params, compute_tw_wr = compute_tw_wr)
    
    return DNA, dyads, nucl
    

def GetNuc(norigin = np.zeros(3), nframe=np.eye(3), wrapped_bp = 0, nuc_type = nuc_name):
    '''
    Obtain the coordinates of a nucleosome from a pdb file

    Parameters
    ----------
    norigin : ndarray of (3)
        center of mass of the nucleosome
    nframe : ndarray of (3,3)
        frame of the nucleosome (x points to dyad, z points along cylinder axis)
    wrapped_bp : int
        controls the amount of nucleosome wrapped DNA
        negative values indicate #bp unwrapped at each side
        positive values indicate total wrapped #bp distributed symetrically 
    Returns
    -------
     params2 : ndarray of (N, 6)
        The basepair step parameters of DNA
     o2 : ndarray of (3)
         start coordinate of DNA
     f2 : ndarray of (3, 3)
         start frame of DNA
     c_dyad : ndarray of (4, 3)
         coordinate frame of the dyad basepair
     
    '''     
    n=nuc.NucPose()
    n.from_file(nuc_type)
    dna = HelixPose(n.params)
    i = 0
    n_bp = dna.n_bp

    o_nuc, f_nuc = nuc.get_nframe(dna.coord, n.d_index, fixed)
    P = nuc.of2coords(o_nuc, f_nuc)
    Q = nuc.of2coords(norigin, nframe)


    transform = nuc.get_transformation(P,Q)
    o_dyad = dna.coord[n.d_index]
    f_dyad = dna.frames[n.d_index]
    c_dyad = nuc.apply_transf_coords(nuc.of2coords(o_dyad, f_dyad) , transform)
        
    c0 = nuc.of2coords(dna.coord[0], dna.frames[0])

    c1 = nuc.apply_transf_coords(c0, transform)
    o1, f1 = nuc.coords2of(c1)
    dna1 = HelixPose(n.params, frame0 = np.transpose(f1))

    
    f2 = dna1.frames[i]
    o2 = o1 + dna1.coord[i]
    params2 = dna1.params[i:i+n_bp]

    return params2, o2, f2, c_dyad
    
def insert_nucs(dna_params, nuc_params, dna_dyads = None,  nuc_dyad = 75, wrapped_bp = None, nuc_type = nuc_name):
    '''
    Fix basepair parameters for nucleosomes in DNA at dyad positions

    Parameters
    ----------
    dna_params : basepair step parameters for DNA
    nuc_params : basepair step parameters for nucleosome
    dna_dyads : dyad positions in DNA (bp)
    nuc_dyad : dyad positions in nucleosome (bp)
    wrapped_bp : controls DNA wrapping around histone octamer(bp)
        negative values indicate #bp unwrapped at each side
        positive values indicate total wrapped #bp distributed symetrically 
    Returns
    -------
     dna_params : ndarray of (N, 6)
        The basepair step parameters of DNA
    free_dna : ndarray of (N)
        indicates wether the basepair is free (1) or constraint in nucleosome (0)     
    '''     

    if dna_dyads is None:
        dna_dyads = [np.array(dna_params.shape[0])/2]
        print dna_dyads

    if wrapped_bp is None or wrapped_bp == 0: bp_index = np.arange(nuc_params.shape[0])-nuc_dyad
    elif wrapped_bp < 0:
        bp_index = (np.arange(nuc_params.shape[0])-nuc_dyad)[-wrapped_bp/2:wrapped_bp/2]

    elif wrapped_bp > 0:
        bp_index = (np.arange(nuc_params.shape[0])-nuc_dyad)[wrapped_bp/2:-wrapped_bp/2]
        
    free_dna = np.ones(dna_params.shape[0])
    for dna_i in dna_dyads:
        for i in bp_index:
            if dna_i + i >= 0 and dna_i + i <= dna_params.shape[0]-1:
                dna_params[dna_i + i,:] = nuc_params[nuc_dyad +i,:]
                free_dna[dna_i + i] = 0
    return dna_params, free_dna


###
###         FUNCTIONS FOR CREATING A CHROMATIN FIBER
###

def set_fiber_params(N, diameter = 330, rise = 100, nld = 17, hand = 'left', ribs = 1):
    """
    Calculate the helical step parameters of a folded chromatin fiber. Step
    parameters are defined relative to nucleosome frames
    
    Parameters
    ----------
    N : number of nucleosomes
    diameter : diameter of the fiber in Angstrom
    rise : distance between nucleosome centers
    nld : nucleosome line density in Angstrom
    hand : chirality of the fiber
    
    Returns
    -------
    fiber_params : nucleosome step parameters of chromatin fiber
    """
    radius = diameter/2.0 - 45.0   
    if ribs == 1:    
        angle = np.pi * np.sqrt(rise**2 - nld**2)/ ( np.pi * radius) #paraxial approximation
        if hand == 'left': angle = -angle
        twist = angle * nld/rise  # Babette
        # twist = angle * rise/nld  # Bert --> wrong
        roll = np.sqrt(angle**2 - twist**2)
        param = np.asarray([0, 0, rise, 0 , roll, twist])
        fiber_params = np.tile(param, (N-1, 1))
    
    if ribs == 2:
 
        cg = 2.0*nld/rise
        theta = np.sqrt(rise**2-4*nld**2)/radius
        znuc1 = np.array([0.0, np.sqrt(1.0-cg**2), cg])
        znuc2 = np.array([np.sqrt(1.0-cg**2)*np.sin(.5*theta), -np.sqrt(1.0-cg**2)*np.cos(.5*theta), cg])          

        
        xnuc1 = np.array([-1.0, 0, 0])
        xnuc2 = np.array([np.cos(.5*theta), np.sin(.5*theta),0])
       
        
        ynuc1 = -np.cross(znuc1, xnuc1)
        ynuc2 = -np.cross(znuc2, xnuc2)

        
        frame1 = np.array([xnuc1, ynuc1, znuc1])
        frame2 = np.array([xnuc2, ynuc2, znuc2])
 
        onuc1 = np.array([radius , 0 , 0])
        onuc2 = np.array([-radius*np.cos(.5*theta), -radius*np.sin(.5*theta), nld])


        
        dr = np.dot(frame1,onuc2-onuc1)
        df = np.dot(frame2,np.transpose(frame1))
        param = util.coords2params(dr,np.transpose(df))

        fiber_params = np.tile(param, (N-1, 1))
    return fiber_params
    
def get_fiber_dyad_frames(fiber_params, nuc_d, nuc_type = nuc_name):
    '''
    Convert the parameters that define a stacked fiber to the
    frame coordinates of the dyad basepairs

    Parameters
    ----------
    fiber_params : ndarray of (N, 6)
        The step paremeters of the stacked fiber

    Returns
    -------
    res : ndarray of (N, 4, 3)
        The coordinates of the dyad frames
    '''  
    fiber = HelixPose(fiber_params)
    res = []
    
    for i in range(len(fiber.frames)):
        n_f = np.transpose(fiber.frames[i])
        n_o = fiber.coord[i]

        dna_params, dna_o, dna_f0, c = GetNuc(norigin = n_o, nframe = n_f, nuc_type = nuc_type)
        n_dna = HelixPose(dna_params, frame0 = dna_f0)

        dyad_coords = nuc.of2coords(n_dna.coord[nuc_d]+dna_o, np.transpose(n_dna.frames[nuc_d]))
        res.append(dyad_coords)
    return res
   
def cast_nucs_fiber(DNA, dyads, dyad_frames, cut = None):
    '''
    fix the location of the nucleosomes in a chromatin fiber
    according the the dyad frames set in dyad_frames. The basepair step 
    parameters are changed in between the nucleosomes to do so. Note that this
    not perfect for large changes at the cut. But as the DNA relaxes, these
    errors will reduce.
    
    The DNA instance will be changed (for time being ;-) )

    Parameters
    ----------
    DNA : instance of HelixMC
    dyads : ndarray of (N)
        bp position of dyads
    dyad_frames : ndarray of (N, 4, 3)
         coordinate frames of the dyad basepairs from stacked helix
    cut : int
        location of the change of step parameters default is halfway 2 nucleosomes
        
    Returns
    -------
     params : ndarray of (N, 6)
        The new basepair step parameters of DNA
     
    '''  
    params = DNA.params    

    for i in range(len(dyads)-1):
# define cut position
        if cut is None:
            cut_bp = dyads[i] + int(0.5*(dyads[i+1]-dyads[i]))
        else:
            cut_bp = dyads[i] + cut
# move end linker nuc[i]
        Q1 = dyad_frames[i]
        P1 = nuc.of2coords(DNA.coord[dyads[i]], np.transpose(DNA.frames[dyads[i]]))
        transform1 = nuc.get_transformation(P1, Q = Q1)
        c_end1 = nuc.of2coords(DNA.coord[cut_bp], np.transpose(DNA.frames[cut_bp]))
        c_end1 = nuc.apply_transf_coords(c_end1, transform1)
        o1, f1 = nuc.coords2of(c_end1)
# move start linker nuc[i+1]
        Q2 = dyad_frames[i+1]
        P2 = nuc.of2coords(DNA.coord[dyads[i+1]],np.transpose(DNA.frames[dyads[i+1]]))
        transform2 = nuc.get_transformation(P2, Q = Q2)
        c_start2 = nuc.of2coords(DNA.coord[cut_bp+1], np.transpose(DNA.frames[cut_bp+1]))
        c_start2 = nuc.apply_transf_coords(c_start2, transform2)   
        o2, f2 = nuc.coords2of(c_start2)
# calc step parameters connection      
        params[cut_bp] = util.frames2params(o1, o2, np.transpose(f1), np.transpose(f2))



    DNA.set_params(params)

    return params
    
def create_curved_linker(n_bp, a):
    '''
    Calculate curved helix parameters as defined by prameters a 

    Parameters
    ----------
    n_bp : int
        linker length
    a : ndarray of (6)
        definition of curved linker
        a[0] = twist (rad)
        a[1] = amplitude
        a[2] = phase (rad)
        a[3] = frequency relative to twist 
        a[4] = modulation frequency
        
    Returns
    -------
     params : ndarray of (N, 6)
        The new basepair step parameters of DNA
     
    '''    
    random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')
    params = np.tile(random_step.params_avg, (n_bp, 1))
    for i in range(n_bp):
        params[i][5] = 2*np.pi/ a[0]   
        params[i][4] =  a[1] * np.cos( a[3]*(2*np.pi/a[0])*float(i) + a[2])
        params[i][4] =  params[i][4] * np.cos(a[4]*2 * np.pi*(float(i)/ n_bp - 0.5))
    return params



##########################
#Calculate Fiber Energies#
##########################

def dna_wrap_unwrap(i, dna, dyad, nucl):
    '''
    Calculate free energy and wrap DNA onto nucleosome if i is a fixed basepair

    Parameters
    ----------
    i : int
        Moving basepair.
    dna : HelixPose
        DNA parameters
    dyad : int
        basepair index of the dyad closest to i
    nucl : NucPose
        parameters of nucleosome template
        
    Returns
    -------
    g_0 : float
        Free energy of wrapping of the nucleosome
    g_1 : float
        Free energy of the nucleosome with nucleosomal DNA rewrapped for i
        up to /from dyad
    params_out : ndarray of (N, 6)
        The new basepair step parameters of DNA with nucleosomal DNA rewrapped
     
    '''
    params_out = dna.params
    new_coords = nuc.get_fixed_coord(dna, dyad, nucl.fixed_i)
    g = nucl.fixed_k*(nucl.fixed_coord-new_coords)**2
    g = np.clip(np.sum(g, axis = 1), 0,  nucl.fixed_g)
    g_0 = np.sum(g)    
    fixed_bp = np.where(nucl.fixed_i + dyad == i)[0]
    
    for n_i in fixed_bp:
        if i <= dyad:
            start = n_i
            end =  np.where(g < nucl.fixed_g)[0][0] + 1
        else:
            start = np.where(g < nucl.fixed_g)[0][-1]
            end = n_i + 1
        if start < end:
            params_out[dyad + nucl.fixed_i[start]:  dyad + nucl.fixed_i[end]+1] = \
               nucl.params[nucl.d_index + nucl.fixed_i[start]-1:  nucl.d_index + nucl.fixed_i[end]]
            g[start : end+1] = 0
    g_1 = np.sum(g)

    return g_0, g_1, params_out
    
    
def dna_wrap_unwrap2(i, dna, dyad, nucl, prev_best):
    '''
    Calc free energy and wrap DNA onto nucleosome if i is a fixed basepair

    Parameters
    ----------
    i : int
        Moving basepair.
    dna : HelixPose
        DNA parameters
    dyad : int
        basepair index of the dyad closest to i
    nucl : NucPose
        parameters of nucleosome template
        
    Returns
    -------
    g_0 : float
        Free energy of wrapping of the nucleosome
    g_1 : float
        Free energy of the nucleosome with nucleosomal DNA rewrapped for i
        up to /from dyad
    params_out : ndarray of (N, 6)
        The new basepair step parameters of DNA with nucleosomal DNA rewrapped
    best_i : int
        the contactpoint that is closed
     '''
    params_out = dna.params
    g_0s = []
    
    for j in range(14):
        new_coords = nuc.get_fixed_coord2(dna, dyad, nucl.fixed_i,fixed_bp = j)
        g = nucl.fixed_k*(nucl.fixed_coord2[j]-new_coords)**2    
        g = np.clip(np.sum(g, axis = 1), 0,  nucl.fixed_g)
        g_0 = np.sum(g)
        g_0s.append(g_0)
        
    g_0s = np.array(g_0s)
    best_is = np.where(g_0s == g_0s.min())    
    
    best_i = np.argmin(g_0s)
    best_i = int(np.mean(best_is))
    
    if len(best_is) == 1:
        best_i = prev_best
    
    new_coords = nuc.get_fixed_coord2(dna, dyad, nucl.fixed_i,fixed_bp = best_i)
    g = nucl.fixed_k*(nucl.fixed_coord2[best_i]-new_coords)**2 
    g = np.clip(np.sum(g, axis = 1), 0,  nucl.fixed_g)
    g_0 = np.sum(g)
    fixed_bp = np.where(nucl.fixed_i + dyad == i)[0] 
    
    for n_i in fixed_bp:
        if i <= nucl.fixed_i[best_i]+dyad:
            start = n_i
            end =  np.where(g < nucl.fixed_g)[0][0]
        else:
            start = np.where(g < nucl.fixed_g)[0][-1]
            end = n_i
        if start < end:
            params_out[dyad + nucl.fixed_i[start]:  dyad + nucl.fixed_i[end]+2] = \
               nucl.params[nucl.d_index + nucl.fixed_i[start]-1:  nucl.d_index + nucl.fixed_i[end]+1]
            g[start : end+1] = 0
        
    g_1 = np.sum(g)

    return g_0, g_1, params_out, best_i
    
    
def dna_type(i, dna, dyads, nucl):
    '''
    Determine the constraints on basepair i

    Parameters
    ----------
    i : int
        Moving basepair.
    dna : HelixPose
        DNA parameters
    dyads : ndarray of (N)
        basepair index of the dyads
    nucl : NucPose
        parameters of nucleosome template
        
    Returns
    -------
    dna_type : ndarray of string
        Indicator of constrains of dna basepair
    closest_dyad : int
        The dyad that is closest to the basepair
     
    '''
    dna_type = 'free'
    dyad_index = -1
    for dyad in dyads:
        if i >= dyad+nucl.fixed_i[0] and i < dyad+nucl.fixed_i[-1]+1:
            dna_type = 'nucleosomal'
            dyad_index = np.where(dyads == dyad)[0][0]
    if len(dyads) > 1 and dna_type == 'free':
        if i >= dyads[0] and i < dyads[-1]:
            dna_type = 'linker'
            dyad_index = np.where(dyads > i)[0][0]-1
    closest_dyad = dyads[dyad_index]

    return dna_type, closest_dyad



def MC_move(i, scorefxn, dna, dyads, nucl, best_is, zigzag = False):
    '''
    Workhorse function for Monte Carlo simulation. This version takes into
    account the wrapping/unwrapping of nucleosomal DNA

    Parameters
    ----------
    i : int
        Moving basepair.
    scorefxn : instance of Cscore
    dna : HelixPose
        DNA that is being modified
    dyads : ndarray (N)
        basepair indeces of the dyad closest to i
    nucl : NucPose
        parameters of nucleosome template
    fiberparams : array of size (N,6)
        step parameters of chromatin fiber
    best_is : dictionary
        best fixed basepairs per nucleosome (dyad number is key)
    zigzag : bool, optional
        True if the fiber is a zigzag type
        
    Returns
    -------
    -
     
    '''

    dna_const, dyad_i = dna_type(i, dna, dyads, nucl)
    prev_best = best_is[dyad_i]
    prev_score = scorefxn(dna)
    

    
        

    prev_params = dna.params
    if dna_const == 'nucleosomal':  
        g_0, g_rewrap, params_out, prev_best = dna_wrap_unwrap2(i, dna, dyad_i, nucl, prev_best)
        prev_score = prev_score + g_0

    if len(dyads) > 1:
        NRL = dyads[1]-dyads[0]
        if zigzag:
            dyad_h = None
            dyad_k = None
            if i < dyads[-1] and i >= dyads[0]:
                if dyad_i > i:
                    dyad_j = dyad_i - NRL                   
                    if dyad_i < dyads[-1]:
                        dyad_h = dyad_i + NRL
                    if dyad_j > dyads[0]:
                        dyad_k = dyad_j - NRL
                    fiberparams = MakeZigzagFiberHelix(dna, [dyad_k, dyad_j, dyad_i, dyad_h], best_is)
                
                
                else:
                    dyad_j = dyad_i + NRL
                    if dyad_i > dyads[0]:
                        dyad_h = dyad_i - NRL
                    if dyad_j < dyads[-1]:
                        dyad_k = dyad_j + NRL
                    fiberparams = MakeZigzagFiberHelix(dna, [dyad_h, dyad_i, dyad_j, dyad_k], best_is)
                prev_score = prev_score + FiberEnergy(fiberparams[0])+FiberEnergy(fiberparams[1])


        else:
            if i < dyads[-1] and i >= dyads[0]:
                if dyad_i > i:
                    dyad_j = dyad_i - NRL
                    fiberparams = MakeFiberHelix(dna, [dyad_j,dyad_i], best_is)
                    
                else:
                    dyad_j = dyad_i + NRL
                    fiberparams = MakeFiberHelix(dna, [dyad_i,dyad_j], best_is)
#                print nuci

            
                prev_score = prev_score + FiberEnergy(fiberparams[0])
        
    next_step_data = random_step()    
    dna.update_trial(i, *next_step_data)
    dna.accept_update()
    score = scorefxn(dna)
    
    
        

    if dna_const == 'nucleosomal':
# compare with rewrapped nucleosome   
        no_rewrap_params = dna.params
        
        g_0, g_rewrap, params_out, new_best = dna_wrap_unwrap2(i, dna, dyad_i, nucl, prev_best)
        dna.set_params(params_out)
        rewrap_score = scorefxn(dna)
        delta_g= (rewrap_score + g_rewrap) - (score + g_0)      
        rnd = random.random_sample()
        accept = (delta_g <= 1e-9 or math.exp(-delta_g / kBT) >= rnd)
        if accept:
            score = rewrap_score + g_rewrap
            best_is[dyad_i] = new_best
        else:
            dna.set_params(no_rewrap_params)
            score = score + g_0

    if len(dyads) > 1:
        if zigzag:
            if i < dyads[-1] and i >= dyads[0]:
                if dyad_i > i:
                    fiberparams = MakeZigzagFiberHelix(dna, [dyad_k, dyad_j, dyad_i, dyad_h], best_is)
                else:
                    fiberparams = MakeZigzagFiberHelix(dna, [dyad_h, dyad_i, dyad_j, dyad_k], best_is)
                score = score + FiberEnergy(fiberparams[0])+FiberEnergy(fiberparams[1])
                    
        else:
            if i < dyads[-1] and i >= dyads[0]:
                if dyad_j > dyad_i:
                    fiberparams = MakeFiberHelix(dna, [dyad_i,dyad_j], best_is)
                else:
                    fiberparams = MakeFiberHelix(dna, [dyad_j,dyad_i], best_is)
    
                score = score + FiberEnergy(fiberparams[0])

# compare with previous score
    delta_g = score - prev_score

    rnd = random.random_sample()
    accept = (delta_g <= 1e-9 or math.exp(-delta_g / kBT) >= rnd)
    if not accept:
        dna.set_params(prev_params)
    
    acc = 0
    if accept: acc = 1
    return acc

def MC_move_only_link(i, scorefxn, dna):
    '''
    Workhorse function for Monte Carlo simulation. This version takes into
    account the wrapping/unwrapping of nucleosomal DNA

    Parameters
    ----------
    i : int
        Moving basepair.
    scorefxn : instance of Cscore
    dna : HelixPose
        DNA that is being modified
    dyads : ndarray (N)
        basepair indeces of the dyad closest to i
    nucl : NucPose
        parameters of nucleosome template
    fiberparams : array of size (N,6)
        step parameters of chromatin fiber
    best_is : dictionary
        best fixed basepairs per nucleosome (dyad number is key)
    zigzag : bool, optional
        True if the fiber is a zigzag type
        
    Returns
    -------
    -
     
    '''


    prev_score = scorefxn(dna)
    prev_params = dna.params
    next_step_data = random_step()    
    dna.update_trial(i, *next_step_data)
    dna.accept_update()
    score = scorefxn(dna)

# compare with previous score
    delta_g = score - prev_score

    rnd = random.random_sample()
    accept = (delta_g <= 1e-9 or math.exp(-delta_g / kBT) >= rnd)
    if not accept:
        dna.set_params(prev_params)
    
    acc = 0
    if accept: acc = 1
    return acc

def count_unwrap(dna, dyads, nucl):
    '''
    Count the number of unwrapped unwrapped basepairs   (Bert 22/01/2015)

    Parameters
    ----------
    dna : HelixPose
        DNA parameters
    dyad : int array
        basepair indices of the dyads
    nucl : NucPose
        parameters of nucleosome template
        
    Returns
    -------
    n_connected : int
        dictionary of connected basepairs at each nucleosome. Nucleome number
        is given by location of its dyad at the basepair
     
    '''
    n_connected = {}
    for dyad in dyads:
        new_coords = nuc.get_fixed_coord(dna, dyad, nucl.fixed_i)
        g = nucl.fixed_k*(nucl.fixed_coord-new_coords)**2
        g = np.clip(np.sum(g, axis = 1), 0,  nucl.fixed_g)
        n = len(np.where(g>(nucl.fixed_g-0.1))[0])

        n_connected[dyad] = n

    return n_connected
    
def count_unwrap2(dna, dyads, nucl, prev_i):
    '''
    Count the number of wrapped fixed basepairs in a mononucleosome  (Bert 22/01/2015)

    Parameters
    ----------
    dna : HelixPose
        DNA parameters
    dyad : int array
        basepair indices of the dyads
    nucl : NucPose
        parameters of nucleosome template
    prev_i : int
        index of previously fixed basepair
        
    Returns
    -------
    n: int
        number of connected basepairs at the first nucleosome
     
    '''
    dyad = dyads[0]
    g_0s = []
    for j in range(14):
        new_coords = nuc.get_fixed_coord2(dna, dyad, nucl.fixed_i,fixed_bp = j)
        g = nucl.fixed_k*(nucl.fixed_coord2[j]-new_coords)**2    
        g = np.clip(np.sum(g, axis = 1), 0,  nucl.fixed_g)
        g_0 = np.sum(g)
        g_0s.append(g_0)
        
    g_0s = np.array(g_0s)
    best_i = np.argmin(g_0s)
    if best_i == 0:
        best_i = prev_i
 
 
    new_coords = nuc.get_fixed_coord2(dna, dyad, nucl.fixed_i,fixed_bp = best_i)
    g = np.sum(nucl.fixed_k*(nucl.fixed_coord2[best_i]-new_coords)**2, axis = 1) 
   
    n = g<nucl.fixed_g
    n = n.astype(int)    
    return n


def get_step_list(nucl):
    '''
    Get the list of nucleotide-nucleotide steps of the nucleosome for sequence 
    dependent modeling
    
    Parameters
    ----------
    nucl : instance of NucPose
    
    Returns
    -------
    step_list : list
    '''
    seq = nucl.chains['DNA'][1]
    step_list = []
    for i in range(len(seq)-1):
        step_list.append(seq[i]+seq[i+1])
    return step_list
        
def Boltzmann_Energy(bp_step, step_type = None):
    # type: (object, object) -> object
    '''
    Get the bending energy of a dna basepair using the gaussian approximation
    
    Parameters
    ----------
    bp_step : ndarray of 6
        the basepair step parameters
    step_type : string, optional
        the type of nucleotide-nucleotide interaction
    
    Returns
    -------
    energy : float
        bending energy of the basepair
    '''
    
    if step_type == None:
        DNApars = "DNA_gau.npy"
        file_working = util.locate_data_file(DNApars)
        params = np.load(file_working)
        avg = params[0]
        cov = params[1:]
        prec = np.linalg.inv(cov)
        return 0.5*kBT*np.dot((bp_step-avg),np.dot(prec,(bp_step-avg)))
    
    else:
        key = step_type
        DNApars = "DNA_2.0_noprot.npz"
        file_working = util.locate_data_file(DNApars)
        params = np.load(file_working)[key]
        avg = np.mean(params, axis = 0)
        cov = np.cov(np.transpose(params))
        prec = np.linalg.inv(cov)
        return 0.5*kBT*np.dot((bp_step-avg),np.dot(prec,(bp_step-avg)))
        
        
def CreateFiber(DNA, dyad1, dyad2, n = 5):
    '''
    Create a fiber from a dinucleosome by copying the nucleosome
    Parameters
    ----------
    DNA : instance of HelixPose
    dyad1 : int
        bp index of first dyad
    dyad2: int
        bp index of second dyad
    n : int
        number of nucleosomes in the  final fiber
        
    Returns
    -------
    '''
    pars = DNA.params
    newpars = pars[:dyad1]
    repeat = pars[dyad1:dyad2]
    end = pars[dyad2:]
    NRL = dyad2-dyad1

    dyads = np.arange(dyad1, dyad1+(n)*NRL,NRL)
    for i in xrange(n-1):
        newpars = np.concatenate((newpars,repeat))
    newpars = np.concatenate((newpars,end))
    DNA.set_params(newpars)
    return dyads
    
def MakeFiberHelix(dna, dyads, best_is, nuc_type = nuc_name):
    '''
    Get step parameters for the solenoid chromatin fiber from the DNA pose
    
    Parameters
    ----------
    dna : HelixPose
        the Pose of the DNA
    dyads : array of integers size (N)
        bp indices of the dyads
    best_is : dictionary
        dictionary with the fixed basepair connected to the histones (indices are
        dyad indices)
    Returns
    -------
    Fiberparams : array of floats size (N,6)
        the helix parameters of the chromatin fiber
    '''

    fiberparams = []
    n = nuc.NucPose()
    n.from_file(nuc_type)
    rprev = np.array([0,0,0])
    fprev = np.eye(3)
    NRL = dyads[1]-dyads[0]
    
    for dyad in dyads:
        index = dyad + n.fixed_i[best_is[dyad]]
        n_origin = dna.coord[index] + np.dot(dna.frames[index], n.tvs[best_is[dyad]])
        n_frame = np.dot(dna.frames[index],n.mats[best_is[dyad]])
        dr = np.dot(np.transpose(fprev), n_origin-rprev)
        df = np.dot(np.transpose(fprev),n_frame)
        rprev = n_origin
        fprev = n_frame        
        fiberparams.append(util.coords2params(dr,df))

    fiberparams.pop(0)
    return np.array(fiberparams)

def MakeZigzagFiberHelix(dna, dyads, best_is):
    '''
    Get step parameters for the zigzag chromatin fiber from the DNA pose
    
    Parameters
    ----------
    dna : HelixPose
    
    dyads : ndarray of integers
    
    Returns
    -------
    fps1 : ndarray of floats
        fiber step parameters of dyads 1, 3, 5 ...
    fps2 : ndarray of floats
        fiber step parameters of dyads 2, 4, 6 ...
    '''
    if dyads[0] == None:
        fps1 = None
        fps2 = MakeFiberHelix(dna, dyads[1::2], best_is)
    if dyads[-1] == None:
        fps1 = MakeFiberHelix(dna, dyads[::2], best_is)
        fps2 = None
    
    if dyads[0] != None and dyads[-1] != None:
        
        d1 = dyads[::2]
        d2 = dyads[1::2]
        fps1 = MakeFiberHelix(dna, d1, best_is)
        fps2 = MakeFiberHelix(dna, d2, best_is)
    return [fps1, fps2]    
    
    
    
def FiberEnergy(fiberparams, av = np.array([0,0,100,0,0.809246,-0.139579])):
    '''
    Calculate the energy of one fiber step
    
    Parameters
    ----------
    fiberparams : array of size (6)
        fiber step parameters
    av : array of size (6)
        lowest energy fiber step parameters
    
    Returns
    -------
    E : float
        the energy of this fiber step in pN*nm
    '''


    if fiberparams is None:
        return 0
    fiberparams = fiberparams[0]
    kxy = 100000000
    kpos = 100*kBT/200.0
    krot = 0
    
    A = np.zeros((6,6))
    for i in range(0,2):
        A[i,i] = kxy
    for i in [3]:
        A[i,i] = kpos
    for i in range(3,6):
        A[i,i] = krot    
    E = 0.5*kBT*np.dot((fiberparams-av),np.dot(A,(fiberparams-av)))
#    if E > 20*kBT:
#        E = 20*kBT
    return E


    
def create_fiber(l_handle, n_nuc, NRL, dna_file = 'DinucleosomePoses\\4qlcNRL167.npz', CL = True):
    '''    
    Create a chromatin fiber with prespecified handle lengths from a dinucleosome
    DNA pose file
    
    Parameters
    ----------
    l_handle : int
        lenght of DNA handle in basepairs
    n_nuc : int
        number of nucleosomes in the fiber
    NRL : int
        nucleosome repeat length in basepairs
    dna_file : string
        name of the saved HelixPose file
        
    Returns
    -------
    dna : instance of HelixPose
    
    dyads : list

    nucl : instance of NucPose
    '''
    dna = HelixPose(params = np.load(dna_file)['params'] , frame0 = np.load(dna_file)['frame0'])
    n_bp = len(dna.coord)
    dna2, dyads, nucl = create_dna(n_bp, 2, NRL, unwrap = -20, nuc_file = (dna_file.split('\\')[1]).split('N')[0]+'.3DNA')
    
    if CL == True:
        n_bp_new = NRL + 147
        diff = n_bp - n_bp_new
        newpars = dna.params[diff/2-1:-(diff/2)]
        for i in range(len(dyads)):
            dyads[i] = nucl.d_index + NRL * i -1
        dna.set_params(newpars)
    dyads = CreateFiber(dna, dyads[0],dyads[1], n_nuc)
#    for i in range(len(dyads)):
#        dyads[i] = dyads[i] + l_handle
    dyads = dyads + l_handle    
    handles = np.tile(random_step.params_avg, (l_handle, 1))
    newpars = dna.params
    
    newpars = np.concatenate((handles,newpars,handles))
    dna.set_params(newpars)
    
    
    
    return dna, dyads, nucl
    
def create_fiber2(l_handle, n_nuc, NRL, dna_file = 'DinucleosomePoses\\4qlcNRL167.npz', nuc_file = '1KX5.3DNA', CL = True):
    '''    
    Create a chromatin fiber with prespecified handle lengths from a dinucleosome
    DNA pose file
    
    Parameters
    ----------
    l_handle : int
        lenght of DNA handle in basepairs
    n_nuc : int
        number of nucleosomes in the fiber
    NRL : int
        nucleosome repeat length in basepairs
    dna_file : string
        name of the saved HelixPose file
        
    Returns
    -------
    dna : instance of HelixPose
    
    dyads : list

    nucl : instance of NucPose
    '''
    dna = HelixPose(params = np.load(dna_file)['params'] , frame0 = np.load(dna_file)['frame0'])
    n_bp = len(dna.coord)
    dna2, dyads, nucl = create_dna(n_bp, 2, NRL, unwrap = -20, nuc_file = nuc_file)
    
    if CL == True:
        n_bp_new = NRL + 147
        diff = n_bp - n_bp_new
        newpars = dna.params[diff/2-1:-(diff/2)]
        for i in range(len(dyads)):
            dyads[i] = nucl.d_index + NRL * i -1
        dna.set_params(newpars)
    dyads = CreateFiber(dna, dyads[0],dyads[1], n_nuc)
#    for i in range(len(dyads)):
#        dyads[i] = dyads[i] + l_handle
    dyads = dyads + l_handle    
    handles = np.tile(random_step.params_avg, (l_handle, 1))
    newpars = dna.params
    
    newpars = np.concatenate((handles,newpars,handles))
    dna.set_params(newpars)
    
    
    
    return dna, dyads, nucl
    
def plotdna(DNA, dyads, plt_range = 50, nuc_type = nuc_name):
    '''
    Plot DNA and show the location of the dyads

    Parameters
    ----------
    DNA : instance of HelixPose
        the DNA to be plotted
    dyads : list
        list of bp indices of the dyads
    plt_range : int, optional
        xyz range (nm), default is 50 nm

      '''

    #plt.close()
    rb_width=2.0
    rb_width=2.0
    color='kb' 
    
    rb_vec = DNA.rb_vec
    coord = DNA.coord/10
    n_centers = []
    vecss = []
    num = 10
    nuc_dum = nuc.NucPose()
    nuc_dum.from_file(nuc_type)
    for d in dyads:
        index = nuc_dum.d_index+fixed[num]
        d1 = (nuc_dum.n_coords[0]-nuc_dum.coords[index])/10
        d2 = np.dot(np.transpose(nuc_dum.frames[index]),d1)
        n_center = coord[d+fixed[num]] + np.dot(DNA.frames[d+fixed[num]],d2)    
        n_centers.append(n_center)
    
    
        n_frame = np.dot(np.dot(DNA.frames[d+fixed[num]],np.transpose(nuc_dum.frames[index])),np.transpose(nuc_dum.n_frame))
    
        vecs =  [n_frame[:,0], n_frame[:,1], n_frame[:,2]]
        vecss.append(vecs)
    kleur = 'bgr'
    
    
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    for j in range(len(dyads)):
        for i in range(3):
            n_center = n_centers[j]
            vecs = vecss[j]
            ax.scatter3D(n_center[0], n_center[1], n_center[2])
            ax.plot3D([n_center[0], n_center[0] + 10*vecs[i][0]],\
                      [n_center[1], n_center[1] + 10*vecs[i][1]],\
                        [n_center[2], n_center[2] + 10*vecs[i][2]], kleur[i])
    
    ax.scatter3D(coord[dyads,0],coord[dyads,1],coord[dyads,2])
    
    ax.plot3D(coord[:, 0], coord[:, 1], coord[:, 2], color[0] + '-')


    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')
    
    if plt_range is None:
        plt_range = DNA.n_bp*0.34      
    ax.set_xlim3d(-plt_range/2,plt_range/2)
    ax.set_ylim3d(-plt_range/2,plt_range/2)
    ax.set_zlim3d(0,plt_range)


def plotdna2(DNA, dyads, nucl, plt_range = 500):
    '''
    Plot DNA and show the location of the dyads

    Parameters
    ----------
    DNA : instance of HelixPose
        the DNA to be plotted
    dyads : list
        list of bp indices of the dyads
    plt_range : int, optional
        xyz range (nm), default is 50 nm

      '''

    #plt.close()
    rb_width=2.0
    rb_width=2.0
    color='kb' 
    
    rb_vec = DNA.rb_vec
    coord = DNA.coord
    NRL = 198
    #NRL = dyads[1]-dyads[0]
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.scatter3D(coord[dyads,0],coord[dyads,1],coord[dyads,2])
    ax.plot3D(coord[:, 0], coord[:, 1], coord[:, 2], color[0] + '-')
    color2 = 'kbgr'
    a = 100
    for d in dyads:
        fixed_frame = np.transpose(DNA.frames[d+nucl.fixed_i[7]-NRL%2])
        fixed_origin = coord[d+nucl.fixed_i[7]-NRL%2]
        Q = nuc.of2coords(fixed_origin, fixed_frame)
        P = nuc.of2coords(nucl.coords[nucl.d_index+nucl.fixed_i[7]], \
                            np.transpose(nucl.frames[nucl.d_index+nucl.fixed_i[7]]))
        tf = nuc.get_transformation(P,Q)
        n_coord = nuc.apply_transf_coords(nucl.n_coords,tf)
        for i in range(1,4):
            ax.plot3D([n_coord[0,0],n_coord[0,0]+a*(n_coord[i,0]-n_coord[0,0])],\
            [n_coord[0,1],n_coord[0,1]+a*(n_coord[i,1]-n_coord[0,1])],\
            [n_coord[0,2],n_coord[0,2]+a*(n_coord[i,2]-n_coord[0,2])], color2[i]+'-')
        ax.scatter3D(coord[nucl.fixed_i +d, 0],coord[nucl.fixed_i +d, 1],coord[nucl.fixed_i +d, 2])

    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')

    if plt_range is None:
        plt_range = DNA.n_bp*0.34      


    ax.set_xlim3d(-plt_range/2,plt_range/2)
    ax.set_ylim3d(-plt_range/2,plt_range/2)
    ax.set_zlim3d(0,plt_range)

    plt._show()

# # are we gonna assume this works?
# a = set_fiber_params(2, diameter = 1000.0, rise = 100.0, nld = 100, ribs = 1)
# print(a)
