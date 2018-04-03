# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:43:39 2015

@author: John

edit: bert

Add functionalities for nucleoosmes to HelixMC
"""
import numpy as np
from helixmc.pose import HelixPose
import os, sys
from helixmc import util
from helixmc import kBT
from matplotlib import pyplot as plt


def find(lst, predicate):
    return (i for i, j in enumerate(lst) if predicate(j)).next()
    

def coords_add_offset(coords, offset):
    """
    Add offset to coordinates
    
    Parameters
    ----------
    coords : ndarray of (N,3)
        coordinates
    offset : ndarray of (3)
        origin coordinates

    
    Returns
    -------
    coords_out: new coordinates    

    """
    N = coords.shape[0]
    coords_out = coords - np.tile(offset, (N, 1))
    return coords_out
    

def of2coords(o, f):
    """
    converts origin and frames to coordinates
    
    Parameters
    -------------
    o : ndarray of (3), origin coordinates
    f : ndarray of (3,3), frame coordinates
    
    Returns
    ---------
    c: coordinates of origin and endpoint of frame
    """
    c = [o]
    for i in range(len(f[0])):
        c.append(f[i]+o)
        
    return np.asarray(c)
    
def coords2of(c):
    """
    convert coordinates to origin and frame, c[0] is the origin 
    
    Parameters
    ----------
    c: coordinates of origin and endpoints of frame    

    
    Returns
    -------
    o : ndarray of (3)
        origin coordinates
    f : ndarray of (3,3)
        frame coordinates
    """    
    o = c[0]
    f = []
    for i in range(3):
        f.append(c[i+1]-o)
    return o, np.asarray(f)

def apply_transf_coords(coords, trans_def):
    """
    Apply rigid transformation to coords using transformation
    parameters obtained with Kabsch algorithm
    
    Parameters
    ----------
    coords : ndarray of (N,3)
    trans_def : list containing center of rotation, rotation matrix and
                center after rotation
    
    Returns
    -------
    coords_out: ndarray of (N,3)    
    """
    N = coords.shape[0]
    coords_out = coords - np.tile(trans_def[0], (N, 1))
    coords_out = np.dot(coords_out, trans_def[1])
    coords_out = coords_out + np.tile(trans_def[2], (N, 1))
    
    return np.asarray(coords_out)
    
def apply_transf_of(o1, f1, trans_def):
    """
    Apply rigid transformation to coords using transformation
    parameters obtained with Kabsch algorithm
    
    Parameters
    ----------
    o1 : ndarray of (3)
        origin coordinates
    f1 : ndarray of (3,3)
        frame
    trans_def : list containing center of rotation, rotation matrix and
                center after rotation
    
    Returns
    -------
    [o2, f2] : ndarray of (3)
        origin coordinates and frame
    """
    c1 = of2coords(o1,f1)
    c2 = apply_transf_coords(c1, trans_def)
        
    return coords2of(c2)

def get_transformation(P, Q = None):
    """
    Align de coordinates defined by P such that
    they will overlap with Q, using Kabsch algorithm
    See http://en.wikipedia.org/wiki/Kabsch_algorithm
    
    Parameters
    ----------
    P : ndarray of (N,3)
        current coordinates
    Q : ndarray of (N,3)
        target coordinates
        if None: Q = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
    
    Returns
    -------
    Pc: translation vector before rotation    
    U : Rotation matrix
    Qc : translation vector after rotation

    """
    if Q is None:
        Q = np.asarray([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
    Pc = sum(P)/(1.0*len(P))
    # print Pc
    Qc = sum(Q)/(1.0*len(Q))
    # print Qc
    C = np.dot(np.transpose(P - Pc), Q - Qc)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return Pc, U, Qc
    
    
def get_nframe(coords, dyad, fixed):
    """
    Calculate the center of mass and the reference frame of a nucleosome
    The nucleosome must have at least 40 bp wrapped on both sides of the dyad
    
    Parameters
    ----------
    coords : coordinates nuclesomal DNA in Angstroms
    dyad : central basepair of nucleosome

    Returns
    -------
    origin : ndarray of (3)
    frame : ndarray of (3,3)
    """

    fixed_i = fixed + dyad
    frame = np.zeros((3,3))


    cm = np.mean(coords[fixed_i], axis = 0)

    Nx = coords[dyad]-cm
    Nx = Nx/np.linalg.norm(Nx)

    Nz = np.mean(coords[fixed_i[:7],:], axis = 0)-np.mean(coords[fixed_i[7:],:], axis = 0)
    Nz = Nz/np.linalg.norm(Nz)

    Ny = np.cross(Nx, Nz)
    Ny = Ny/np.linalg.norm(Nz)
    
    Nz = np.cross(Nx, Ny)
    Nz = Nz/np.linalg.norm(Nz)
    origin = cm
    frame = np.array([Nx, Ny, Nz])
    return origin, frame
    

    
def seq3_to_1(seq3):
    '''
    Turn a three letter protein into a one letter protein.
    Workd also for DNA 'DA ' strings

    Parameters
    ----------
    seq3 : str(4) aminoacid or DNA code
    
    Returns
    -------
    seq1 : str

    '''
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X',
     'DA': 'A', 'DC': 'C','DG': 'G', 'DT': 'T', '   ':''}
    seq3 = seq3.strip()
    seq1=''
    for i in range(0, len(seq3), 4):
        seq1 += d[seq3[0+i:3+i].strip()]
    return seq1

    
def ReadPdb(pdb_file):
    '''
    Get Chains, Chain type and coords of proteins from pdb file 

    Parameters
    ----------
    pdb_file : string
    
    Returns
    -------
    chain_dict: {chain, type, coords}

    '''
#    pdb_file = "P:\\My Documents\\MasterProjectBert\\ChromatinMC\\data\\1KX5.pdb"
    pdb_file = "PDBs\\1KX5.pdb"

    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    cpd = [l for l in pdb if l[0:6] == "COMPND"]
    chains=[]
    molecules=[]
    keywords = ['DNA', 'H2A', 'H2B', 'H3', 'H4']
    for l in cpd:
        s =l[11:].split(':')
        if s[0] == 'MOLECULE':
            for k in keywords:
                if s[1].find(k) >= 0: 
                    if k in chains:
                        chains.append(k + '*')
                    else: chains.append(k)
        if s[0] == 'CHAIN':
            s=s[1].split(';')[0].split(',')
            i=0
            for m in s:
                molecules.append(m.lstrip())
                if i > 0: 
                    chains.append(chains[-1]+'*')
                i+=1

    chain_dict = dict([(c,['','', np.zeros((1,3)), np.zeros(1)]) for c in chains])
    i = 0
    for c in chains:
        chain_dict[c][0] = molecules[i]
        i += 1
    chain_dict.pop('DNA*')
    
    # Use SEQRES for sequence 
    seq = [l for l in pdb if l[0:6] == "SEQRES"]
    for l in seq:
        for i in chain_dict:
            if chain_dict[i][0] == l[11]:
                chain_dict[i][1] += l[19:71]

    for i in chain_dict:
        chain_dict[i][1] = seq3_to_1(chain_dict[i][1])
        chain_dict[i][2] = np.zeros([len(chain_dict[i][1]),3])
    #removed /0 from "chain_dict[i][2] = np.zeros([len(chain_dict[i][1]),3])/0"
    #-bert
    
    # Grab all CA from ATOM entries
    # Grab B -factor for Phosphates
    B_factor = np.zeros(len(chain_dict['DNA'][2])-1)
    P_I = np.zeros((len(chain_dict['DNA'][2])-1,3))

    
    start_i = None  
    for l in pdb:
        if l[0:6] == "ATOM  " and l[13:16] == "CA ":
            nc = np.fromstring(l[30:55], sep = ' ')
            for i in chain_dict:
                if chain_dict[i][0] == l[21]:
#                    print int(l[22:26]
                    chain_dict[i][2][int(l[22:26])-1,:] = nc         
        if l[0:6] == "ATOM  " and l[13:16] == "P  ":
            if start_i == None:
                start_i = int(l[22:27])
            B_factor[int(l[22:27])-start_i] += float(l[61:66])
            if l[21]=='I':
                P_I[int(l[22:27])-start_i] = np.fromstring(l[30:55], sep = ' ')

#    chain_dict['DNA'][3]=B_factor
    av_I = np.mean(P_I, axis = 0)

    dI = np.sqrt(np.sum((P_I-av_I)**2, axis = 1))
    chain_dict['DNA'][3]= dI
#    plt.plot(np.arange(len(dI)), dI)
#
#    plt.plot(np.arange(len(B_factor)), B_factor)
#    plt.show()
    return chain_dict

def crossproduct(u,v):
    w1 = u[1]*v[2]-u[2]*v[1]
    w2 = u[2]*v[0]-u[0]*v[2]
    w3 = u[0]*v[1]-u[1]*v[0]
    return np.array([w1,w2,w3])
    
def get_frame(DNA, i):
    base = DNA.coord[i]
    yn = DNA.rb_vec[i]/np.linalg.norm(DNA.rb_vec[i])
    zn = DNA.dr[i]/np.linalg.norm(DNA.dr[i])
    xn = crossproduct(yn, zn)
    x = np.asarray([base, base + xn, base + yn, base + zn])
    return x


    
def get_fixed_coord(dna, dyad, fixed_i, T = False):
     
    coords = []    
    center = dyad + fixed_i[6]
    for i in fixed_i:
        coords.append(dna.coord[i+dyad])
    P = of2coords(dna.coord[center], np.transpose(dna.frames[center]))
    transform = get_transformation(P)
    fixed_coords = apply_transf_coords(np.asarray(coords), transform)
    return fixed_coords          
 
    
def get_fixed_coord2(dna, dyad, fixed_i, fixed_bp = 6):
    
    coords = []    
    center = dyad + fixed_i[fixed_bp]
    for i in fixed_i:
        coords.append(dna.coord[i+dyad])
    P = of2coords(dna.coord[center], np.transpose(dna.frames[center]))
    transform = get_transformation(P)
    fixed_coords = apply_transf_coords(np.asarray(coords), transform)
    
    '''
    center_origin = nucl.DNAcoords[center+nucl.d_index-dyad]
    center_frame = nucl.nucDNAcopy.frames[center+nucl.d_index-dyad]
    center_coords = of2coords(center_origin, center_frame)
    for i in fixed_i:
        coords.append(dna.coord[i+dyad])
    
    P = of2coords(dna.coord[center], np.transpose(dna.frames[center]))
    transform = get_transformation(P, Q = center_coords)
    fixed_coords = apply_transf_coords(np.asarray(coords), transform)
    '''
    return fixed_coords        
 
    

        
class NucPose(object):
    '''
    NucPose object storing the state info of the nucleosome.

    Parameters
    ----------
    
    Attributes
    ----------
    step_list : list
        List of nucleotide-nucleotide steps (eg AC, AG etc)    
    n_bp : int
        Number of base-pairs in the nucleosome.
    params : ndarray of (N,6)
        List of all bp-step parameters of the DNA in Angstrom. 
        A N bp nucleosome has (N-1)*6 parameters.
    d_index : int
        index of dyad bp
    d_origin : ndarray of (3)
        Origin of the dyad in Angstrom
    d_frame : ndarray of (3,3)
        The frame of dyad
    d_coords : ndarray of (4, 3)
        coordinates of dyad frame
    n_origin : ndarray of (3)
        The center of mass of the nucleosome in Angstrom
    n_frame : ndarray of (3,3)
        The frame of nucleosome
    n_coords : ndarray of (4, 3)
        coordinates of nucleosome frame
    coords : ndarray of (N,3)
        The coordinates of the basepairs
    frame0: ndarray of (3,3)
        start frame of pdb DNA coords
    origin0: ndarray of  (3,1)
        start coordinates of pdb DNA coords
    chains : dictionary
        dictionary of all chains {name, sequence, coords}
    fixed_i : list
        indices of fixed basepairs relative to dyad
    l_coords : ndarray of (4,3)
        The coordinates of Glu61 (H2A) and Asp24 (H4) that mediate nucleosome-
        nucleosome interactions through the tail of H4
    '''
    def __init__(self):
        return

    @classmethod
    def from_file(self, input_file):
        '''
        Load pose data from an input file.

        Parameters
        ----------
        input_file : str
        input file (.txt) obtained from w3DNA.

        '''
        self.nuc_type = input_file.split('.')[0]
        if input_file == None:
            input_file = "3LZ0.3DNA"
            print('Default file: ' + input_file)
#        directory = "C:\\Users\\lion\\Desktop\\ChromatinMC\\PDBs\\"
        directory = "PDBs\\3LZ0.3DNA"

#   read protein coords from pdb file

        chains = ReadPdb(directory+ input_file.replace('3DNA','pdb'))
#   get fixed frames
        B = chains['DNA'][3]

        acids = chains['DNA'][1]
        self.step_list = []
        for i in range(len(acids)-1):
            self.step_list.append(acids[i]+acids[i+1])
      
        fixed = []
        for i in range(14):
            j = np.argmin(B)
            fixed.append(j)
            B[j-6:j+6]=np.max(B)

        fixed = np.sort(fixed)
        self.d_index = int(round(np.mean(fixed)))
        self.fixed_i = fixed - self.d_index

        
#   get link coordinates Glu61 (H2A) and Asp24 (H4)
        int_dict={'H2A':60, 'H2A*':60, 'H4':23, 'H4*':23}
        self.l_coords = []
        for locus in int_dict:
            self.l_coords.append(chains[locus][2][int_dict[locus]])
        self.l_coords = np.asarray(self.l_coords)
        
#   read 3DNA file         
#        with open(directory+input_file) as f:
        with open(directory) as f:
            lines = f.read().splitlines()
            
#   Basepair step parameters            
        i=find(lines, lambda x:
        '    step       Shift     Slide      Rise      Tilt      Roll     Twist'
        in x)+1

        j=find(lines[i:], lambda x: '~~~~~~~~~~~' in x)

        self.params=np.zeros((j,6))
        for k in range (0,j):
            if not('----' in lines[i+k]):
                tmp = np.asarray(map(float,lines[i+k].split()[2:]))
            self.params[k:]= np.append(tmp[0:3], np.radians(tmp[3:6]))
        self.n_bp = len(self.params[:,0])+1
        

        
 
        
            
            
        
        

#   Centers of basepairs
        i=find(lines, lambda x: 
        '      bp        Ox        Oy        Oz        Nx        Ny        Nz'
        in x)+1

        coords = np.zeros((j,3))
        frames = np.zeros((j,3,3))
        for k in range (0,j):
            tmp = np.asarray(map(float,lines[i+k].split()[2:]))
            coords[k,:] = tmp[0:3]
            frames[k][0] = tmp[3:6]/np.linalg.norm(tmp[3:6])
            if k is self.d_index:
                self.d_origin = tmp[0:3]
                self.d_frame = tmp[3:6]/np.linalg.norm(tmp[3:6])
        chains['DNA'][2] = coords
        self.chains = chains
        
        
        if self.n_bp > 147:
            loose_ends = 147 - (self.fixed_i[-1]-self.fixed_i[0]+1)
            if loose_ends%2 == 1:
                l1 = loose_ends/2
                l2 = l1 + 1
            else:
                l1 = loose_ends/2
                l2 = l1
            start = self.fixed_i[0] - l1 + self.d_index
            end = self.fixed_i[-1] + l2 +self.d_index
            self.params = self.params[start:end]
            coords = coords[start:end]
            self.d_index = self.d_index-start
#   Dyad frame coords
        nucDNA = HelixPose(self.params)        
        self.d_coords = of2coords(nucDNA.coord[self.d_index], np.transpose(nucDNA.frames[self.d_index]))
        self.coords = nucDNA.coord
        self.frames = nucDNA.frames
#   Nucleosome frame
        self.n_origin, self.n_frame = get_nframe(self.coords, self.d_index, self.fixed_i)
        self.n_coords = of2coords(self.n_origin, self.n_frame)
        

        
#   origin of nucleosomal DNA
        transform = get_transformation(nucDNA.coord[0:len(coords)], Q = coords)   
        framecoords = of2coords(np.zeros(3), np.eye(3))
        framecoords = apply_transf_coords(framecoords, transform)
        self.origin0, self.frame0 = coords2of(framecoords)

#   fixed coordinates of nucleosome
        self.fixed_coord = get_fixed_coord(nucDNA, self.d_index-1, self.fixed_i)
        self.fixed_coord2 = []
        for i in range(0,14):
            self.fixed_coord2.append(get_fixed_coord2(nucDNA, self.d_index-1, \
            self.fixed_i, i))
            
#   create lists of linear transforms of fixed basepair to center of mass
        self.tvs = []
        self.mats = []
        for i in range(0,14):
            index = self.d_index + self.fixed_i[i]
            #multiply those with the frame of the fixed basepair of the dna to get the 
            #frame and position of the nucleosome center of mass
            self.tvs.append(np.dot(np.transpose(self.frames[index]), self.n_coords[0]-self.coords[index]))
            
            self.mats.append(np.dot(np.transpose(self.frames[index]),np.transpose(self.n_frame)))
        
#   stiffnes of fixed basepairs          
        sd_pos = 0.5 # (A)
        k_pos = kBT / sd_pos**2
        self.fixed_k = k_pos
        self.fixed_g = 10*kBT
        



class NucPose_original(object):
    '''
    NucPose object storing the state info of the nucleosome.

    Parameters
    ----------
    
    Attributes
    ----------
    n_bp : int
        Number of base-pairs in the nucleosome.
    params : ndarray of (N,6)
        List of all bp-step parameters of the DNA in Angstrom. 
        A N bp nucleosome has (N-1)*6 parameters.
    d_index : int
        index of dyad bp
    d_origin : ndarray of (3)
        Origin of the dyad in Angstrom
    d_frame : ndarray of (3,3)
        The frame of dyad
    d_coords : ndarray of (4, 3)
        coordinates of dyad frame
    n_origin : ndarray of (3)
        The center of mass of the nucleosome in Angstrom
    n_frame : ndarray of (3,3)
        The frame of nucleosome
    n_coords : ndarray of (4, 3)
        coordinates of nucleosome frame
    coords : ndarray of (N,3)
        The coordinates of the basepairs
    frame0: start frame of pdb DNA coords
    origin0: start coordinates of pdb DNA coords
    chains : dictionary of all chains {name, sequence, coords}
    fixed_i : indices of fixed basepairs relative to dyad
    l_coords : ndarray of (4,3)
        The coordinates of Glu61 (H2A) and Asp24 (H4) that mediate nucleosome-
        nucleosome interactions through the tail of H4
    '''
    def __init__(self):
        return

    @classmethod
    def from_file(self, input_file):
        '''
        Load pose data from an input file.

        Parameters
        ----------
        input_file : str
        input file (.txt) obtained from w3DNA.

        '''
        
        if input_file == None:
            input_file = "3LZ0.3DNA"
            
            print('Default file: ' + input_file)
        directory = "N:\\Chromatin\\Software\\ChromatinMC\\PDBs"
#   read protein coords from pdb file

        chains = ReadPdb(directory+ input_file.replace('3DNA','pdb'))
#   get fixed frames
        B = chains['DNA'][3]

        acids = chains['DNA'][1]
        self.step_list = []
        for i in range(len(acids)-1):
            self.step_list.append(acids[i]+acids[i+1])
      
        fixed = []
        for i in range(14):
            j = np.argmin(B)
            fixed.append(j)
            B[j-6:j+6]=np.max(B)

        fixed = np.sort(fixed)
        self.d_index = int(round(np.mean(fixed)))
        self.fixed_i = fixed - self.d_index


        
#   get link coordinates Glu61 (H2A) and Asp24 (H4)
        int_dict={'H2A':60, 'H2A*':60, 'H4':23, 'H4*':23}
        self.l_coords = []
        for locus in int_dict:
            self.l_coords.append(chains[locus][2][int_dict[locus]])
        self.l_coords = np.asarray(self.l_coords)
        
#   read 3DNA file         
        with open(directory+input_file) as f:
            lines = f.read().splitlines()
            
#   Basepair step parameters            
        i=find(lines, lambda x: 
        '    step       Shift     Slide      Rise      Tilt      Roll     Twist'
        in x)+1

        j=find(lines[i:], lambda x: '~~~~~~~~~~~' in x)

        self.params=np.zeros((j,6))
        for k in range (0,j):
            if not('----' in lines[i+k]):
                tmp = np.asarray(map(float,lines[i+k].split()[2:]))
            self.params[k:]= np.append(tmp[0:3], np.radians(tmp[3:6]))
        self.n_bp = len(self.params[:,0])+1
        
        
        
#        self.origins, self.frames = util.params2coords(self.params)
#        self.d_origin =self.origins[self.d_index]
#        self.d_frame = self.frames[self.d_index]
        
        
        

#   Centers of basepairs
        i=find(lines, lambda x: 
        '      bp        Ox        Oy        Oz        Nx        Ny        Nz'
        in x)+1

        coords = np.zeros((j,3))
        frames = np.zeros((j,3,3))
        for k in range (0,j):
            tmp = np.asarray(map(float,lines[i+k].split()[2:]))
            coords[k,:] = tmp[0:3]
            frames[k][0] = tmp[3:6]/np.linalg.norm(tmp[3:6])
            if k is self.d_index:
                self.d_origin = tmp[0:3]
                self.d_frame = tmp[3:6]/np.linalg.norm(tmp[3:6])
        chains['DNA'][2] = coords
        self.chains = chains
        
        
#   Dyad frame coords
        nucDNA = HelixPose(self.params)
        self.nucDNAcopy = HelixPose(self.params)
        Qn = coords
        P = nucDNA.coord[0:len(Qn)]
        transform = get_transformation(P, Q = Qn)
        
        
        Pn = of2coords(nucDNA.coord[self.d_index], nucDNA.frames[self.d_index])
        self.d_of = Pn

        self.d_coords = of2coords(nucDNA.coord[self.d_index], np.transpose(nucDNA.frames[self.d_index]))
        self.DNAcoords = nucDNA.coord
        self.DNApose = nucDNA
#   Nucleosome frame
        self.n_origin, self.n_frame = get_nframe(self.DNAcoords, self.d_index, self.fixed_i)
        #self.n_frame = np.transpose(self.n_frame)
        self.n_coords = of2coords(self.n_origin, self.n_frame)
        
        #com frame to dyad_frame transform
        self.com_d_transform = get_transformation(self.n_coords,self.d_coords)
        
#   origin of nucleosomal DNA
        transform = get_transformation(nucDNA.coord[0:len(coords)], Q = coords)   
        framecoords = of2coords(np.zeros(3), np.eye(3))
        framecoords = apply_transf_coords(framecoords, transform)
        self.origin0, self.frame0 = coords2of(framecoords)

#   fixed coordinates of nucleosome
        self.fixed_coord = get_fixed_coord(nucDNA, self.d_index-1, self.fixed_i)
        self.fixed_coord2 = []
        for i in range(0,14):
            self.fixed_coord2.append(get_fixed_coord2(nucDNA, self.d_index-1, \
            self.fixed_i, i))
        
#   stiffnes of fixed basepairs          
        sd_pos = 0.5 # (A)
        k_pos = kBT / sd_pos**2
        self.fixed_k = k_pos
        self.fixed_g = 10*kBT
        
        if self.n_bp > 147:
            loose_ends = self.fixed_i[-1]-self.fixed_i[0]+1
            if loose_ends%2 == 1:
                l1 = loose_ends/2
                l2 = l1 + 1
            else:
                l1 = loose_ends/2
                l2 = l1
            self.DNAcoords = self.DNAcoords[self.d_index-l1+self.fixed_i[0]:self.d_index+self.fixed_i[-1]+l2]

        #return self.params