# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:43:39 2015

@author: John

Add functionalities for nucleoosmes to HelixMC
"""

import matplotlib as mpl

mpl.use(u'TkAgg')
mpl.interactive(False)

import numpy as np
from helixmc.pose import HelixPose
from helixmc.util import frames2params
import matplotlib.pyplot as plt
# ChromatinMC modules:
import FileIO3 as fileio

plt.interactive(False)
np.set_printoptions(formatter={'float': '{: 0.1f}, '.format})

pdb_source_dir = "PDBs\\"
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})


def ofs2params(of1, of2):
    o1 = of1[0]
    f1 = of1[1:, :] - o1
    o2 = of2[0]
    f2 = of2[1:, :] - o2
    f1 = np.transpose(f1)
    f2 = np.transpose(f2)
    params = frames2params(o1, o2, f1, f2)
    return params


def find(lst, predicate):
    return (i for i, j in enumerate(lst) if predicate(j)).next()


def of2axis(of, length=60):
    """
    converts originframe to axis for plotting purposes
    """
    o = of[0]
    f = of[1:] - o
    coords_out = []
    for i in np.arange(0, 3):
        for j in np.linspace(0, (i + 1) * length, length * 1):
            coords_out.append(o + j * f[i])
    return np.asarray(coords_out)


def apply_transformation(coords, trans_def):
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


def get_transformation(start, target=None):
    """
    Align de coordinates defined by P such that
    they will overlap with Q, using Kabsch algorithm
    See http://en.wikipedia.org/wiki/Kabsch_algorithm
    """
    P = start
    Q = target
    if Q is None:
        Q = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    Pc = sum(P) / (1.0 * len(P))
    Qc = sum(Q) / (1.0 * len(Q))
    C = np.dot(np.transpose(P - Pc), Q - Qc)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return Pc, U, Qc


def join_o_f(origin, frame):
    of = []
    of.append(origin)
    for f in np.transpose(frame):
        of.append(origin + f)
    return np.asarray(of)


def get_nuc_of(dna, dyad, nucl):
    """
    Calculate the center of mass and the reference frame (=of) of a nucleosome in dna
    """
    tf = get_transformation(get_of(nucl.dna, nucl.dyad), get_of(dna, dyad))
    n_of = apply_transformation(nucl.of, tf)
    return n_of


def get_of(dna, i):
    of = []
    of.append(dna.coords[i])
    for f in np.transpose(dna.frames[i]):
        of.append(dna.coords[i] + f)
    return np.asarray(of)


def get_wrap_params(dna, dyad, fixed):
    fixed_params = []
    dyad_of = get_of(dna, dyad)
    for i in fixed:
        base_of = get_of(dna, dyad + i)
        fixed_params.append(ofs2params(dyad_of, base_of))
    return np.asarray(fixed_params).reshape((-1, 6))


def seq3_to_1(seq3):
    '''
    Turn a three letter protein into a one letter protein.
    Works also for DNA 'DA ' strings

    Parameters
    ----------
    seq3 : str(4) aminoacid or DNA code
    
    Returns
    -------
    seq1 : str

    '''
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER': '*',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'XAA': 'X',
         'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', '   ': ''}
    seq3 = seq3.strip()
    seq1 = ''
    for i in range(0, len(seq3), 4):
        seq1 += d[seq3[0 + i:3 + i].strip()]
    return seq1


def read_pdb(pdb_file):
    '''
    Get Chains, Chain type and coords of proteins from pdb file 

    Parameters
    ----------
    pdb_file : string
    
    Returns
    -------
    chain_dict: {chain, type, coords}

    '''
    #    print '###>'+ pdb_file
    f = open(pdb_file, 'r')
    pdb = f.readlines()
    f.close()

    cpd = [l for l in pdb if l[0:6] == "COMPND"]
    chains = []
    molecules = []
    keywords = ['DNA', 'H2A', 'H2B', 'H3', 'H4']
    for l in cpd:
        s = l[11:].split(':')
        if s[0] == 'MOLECULE':
            for k in keywords:
                if s[1].find(k) >= 0:
                    if k in chains:
                        chains.append(k + '*')
                    else:
                        chains.append(k)
        if s[0] == 'CHAIN':
            s = s[1].split(';')[0].split(',')
            i = 0
            for m in s:
                molecules.append(m.lstrip())
                if i > 0:
                    chains.append(chains[-1] + '*')
                i += 1

    chain_dict = dict([(c, ['', '', np.zeros((1, 3)), np.zeros(1)]) for c in chains])
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
        chain_dict[i][2] = np.zeros([len(chain_dict[i][1]), 3])
    # removed /0 from "chain_dict[i][2] = np.zeros([len(chain_dict[i][1]),3])/0"
    # -bert

    # Grab all CA from ATOM entries
    # Grab B -factor for Phosphates
    B_factor = np.zeros(len(chain_dict['DNA'][2]) - 1)
    P_I = np.zeros((len(chain_dict['DNA'][2]) - 1, 3))

    start_i = None
    for l in pdb:
        if l[0:6] == "ATOM  " and l[13:16] == "CA ":
            nc = np.fromstring(l[30:55], sep=' ')
            for i in chain_dict:
                if chain_dict[i][0] == l[21]:
                    chain_dict[i][2][int(l[22:26]) - 1, :] = nc
        if l[0:6] == "ATOM  " and l[13:16] == "P  ":
            if start_i == None:
                start_i = int(l[22:27])
            B_factor[int(l[22:27]) - start_i] += float(l[61:66])
            if l[21] == 'I':
                P_I[int(l[22:27]) - start_i] = np.fromstring(l[30:55], sep=' ')

    #    chain_dict['DNA'][3]=B_factor
    av_I = np.mean(P_I, axis=0)

    dI = np.sqrt(np.sum((P_I - av_I) ** 2, axis=1))
    chain_dict['DNA'][3] = dI
    #    plt.plot(np.arange(len(dI)), dI)
    #
    #    plt.plot(np.arange(len(B_factor)), B_factor)
    #    plt.show()
    return chain_dict


def crossproduct(u, v):
    w1 = u[1] * v[2] - u[2] * v[1]
    w2 = u[2] * v[0] - u[0] * v[2]
    w3 = u[0] * v[1] - u[1] * v[0]
    return np.array([w1, w2, w3])


class NucPose(object):
    '''
    NucPose object storing the state info of the nucleosome.

    Parameters
    ----------
    
    Attributes
    ----------
    nuc_type : string
        Name of the pdb file that contains teh coordinates
    step_list : list
        List of nucleotide-nucleotide steps (eg AC, AG etc)    
    dna : HelixMC pose
        nucleosomal DNA in HelixMC framework
    dyad : int
        index of dyad bp
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
        #   read protein coords from pdb file
        filename = pdb_source_dir + input_file
        filename = fileio.change_extension(filename, 'pdb')
        chains = read_pdb(filename)

        #   read 3DNA file
        filename = fileio.change_extension(filename, '3DNA')
        with open(filename) as f:
            lines = f.read().splitlines()

        #   Basepair step parameters
        i = find(lines, lambda x:
        '    step       Shift     Slide      Rise      Tilt      Roll     Twist'
        in x) + 1

        j = find(lines[i:], lambda x: '~~~~~~~~~~~' in x)

        params = np.zeros((j, 6))
        for k in range(0, j):
            if not ('----' in lines[i + k]):
                tmp = np.asarray(map(float, lines[i + k].split()[2:]))
            params[k:] = np.append(tmp[0:3], np.radians(tmp[3:6]))

        #   get fixed frames
        B_factor = chains['DNA'][3]
        fixed = []
        for i in range(14):
            j = np.argmin(B_factor)
            fixed.append(j)
            B_factor[j - 6:j + 6] = np.max(B_factor)

        fixed = np.sort(fixed)
        self.dyad = int(round(np.mean(fixed)))
        self.fixed = fixed - self.dyad

        #   Centers of basepairs
        i = find(lines, lambda x:
        '      bp        Ox        Oy        Oz        Nx        Ny        Nz'
        in x) + 1

        coords = np.zeros((j, 3))
        frames = np.zeros((j, 3, 3))
        for k in range(0, j):
            tmp = np.asarray(map(float, lines[i + k].split()[2:]))
            coords[k, :] = tmp[0:3]
            frames[k][0] = tmp[3:6] / np.linalg.norm(tmp[3:6])
        chains['DNA'][2] = np.asarray(coords)

        # rotate all coords, such that frame0 = origin
        self.dna = HelixPose(params)
        tm = get_transformation(chains['DNA'][2][0:100], self.dna.coords[0:100])
        for chain in chains:
            chains[chain][2] = apply_transformation(chains[chain][2], tm)
        self.chains = chains

        # remove non-nucleosomal DNA
        if len(self.dna.coords) > 147:
            loose_ends = 147 - (self.fixed_i[-1] - self.fixed_i[0] + 1)
            if loose_ends % 2 == 1:
                l1 = loose_ends / 2
                l2 = l1 + 1
            else:
                l1 = loose_ends / 2
                l2 = l1
            start = self.fixed_i[0] - l1 + self.d_index
            end = self.fixed_i[-1] + l2 + self.d_index
            self.d_index = self.d_index - start
            params = params[start:end]
            self.dna = HelixPose(params)

        # get origin and frame of nucleosome
        cm = np.mean(self.dna.coords[self.fixed], axis=0)
        Nx = self.dna.coords[self.dyad] - cm
        Nx = Nx / np.linalg.norm(Nx)
        Nz = np.mean(self.dna.coords[self.fixed[:7], :], axis=0) - np.mean(self.dna.coords[self.fixed[7:], :], axis=0)
        Nz = Nz / np.linalg.norm(Nz)

        Ny = np.cross(Nx, Nz)
        Ny = Ny / np.linalg.norm(Nz)

        Nz = np.cross(Nx, Ny)
        Nz = Nz / np.linalg.norm(Nz)
        origin = cm
        frame = np.array([Nx, Ny, Nz])
        self.of = join_o_f(origin, np.transpose(frame))

        #   get link coordinates Glu61 (H2A) and Asp24 (H4)
        # int_dict = {'H2A': 60, 'H2A*': 60, 'H4': 23, 'H4*': 23}
        # self.l_coords = []
        # for locus in int_dict:
        #     self.l_coords.append(chains[locus][2][int_dict[locus]])
        # self.l_coords = np.asarray(self.l_coords)


def main():
    nuc = NucPose()
    nuc.from_file('1KX5.3DNA')

    coords = []
    for chain in nuc.chains:
        coords.append(nuc.chains[chain][2])

    tf = get_transformation(nuc.of)

    n_coords = []
    for c in coords:
        n_coords.append(apply_transformation(c, tf))

    nuc_ax = apply_transformation(nuc.of, tf)
    n_coords.append(of2axis(nuc_ax, length=60))

    filename = fileio.get_filename(root='1nuc', ext='pov', incr=True)
    fileio.pov_dna(filename, n_coords, range_A=[300, 300], offset_A=[0, 0, 60], show=True, width_pix=500)


if __name__ == '__main__':
    main()
