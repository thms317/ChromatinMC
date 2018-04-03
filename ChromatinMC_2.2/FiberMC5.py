# -*- coding: utf-8 -*-

import matplotlib as mpl

mpl.use(u'TkAgg')
mpl.interactive(False)

import sys
import numpy as np
from helixmc.pose import HelixPose
from helixmc.random_step import RandomStepSimple
from lmfit import Minimizer, Parameters, report_fit, minimize
import matplotlib.pyplot as plt
# ChromatinMC modules:
import FileIO3 as fileio
import NucleosomeMC3 as nMC


plt.interactive(False)
np.set_printoptions(formatter={'float': '{: 0.1f}, '.format})


def create_curved_linker(p=None):
    if p is None:
        p = Parameters()
        p.add('L_bp', value=400)
        p.add('n_nuc', value=2)
        p.add('nld_A', value=100)
        p.add('NRL', value=197)
        p.add('Unwrapped_bp', value=30)
        p.add('Pitch_bp', value=10.380450417509)
        p.add('G_kT', value=71.8913935387115)
        p.add('Phase_rad', value=3.75743188237752)
        p.add('Period', value=0.787822017675862)
        p.add('RelFreq', value=1.09559893586307)

    random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')

    linkerlength_bp = p['NRL'] - 147 + p['Unwrapped_bp']

    params = np.tile(random_step.params_avg, (linkerlength_bp, 1))

    i = np.linspace(-linkerlength_bp // 2, linkerlength_bp // 2 - 1, linkerlength_bp)
    twist = 0 * i + 2 * np.pi / p['Pitch_bp']

    roll = p['Phase_rad'] + p['RelFreq'] * 2 * np.pi * i / p['Pitch_bp']
    roll = np.cos(roll)

    modulation = np.cos(p['Period'] * 2.0 * np.pi * i / linkerlength_bp)
    roll = roll * modulation

    energy = np.sqrt(np.sum(roll ** 2)) / 0.115  # 0.115 = approximate twist stiffness
    # Chou2014PloSCmpBiol: sigma_roll = 5' so k = kT/sigma^2 =4.1/(5*pi/180)^2 = 0.01857
    roll = roll * np.sqrt(p['G_kT']) / energy

    params[:, 5] = twist
    params[:, 4] = roll

    return params


def replace_linker(dna, dyads, linker_params):
    '''
    Replace basepair parameters for linker DNA

    Parameters
    ----------
    dna : basepair step parameters for DNA
    dyads : dyad index in DNA (bp)
    free_bps : array of bps; 0 = nucleosomal DNA, 1 = free DNA
    linker_params : basepair step parameters for linker DNA
    Returns
    -------
    dna_params : ndarray of (N, 6)
        The basepair step parameters of DNA
    free_bps : array of bps; 0 = nucleosomal DNA, 1 = free DNA
    '''
    params = dna.params
    for d1, d2 in zip(dyads, dyads[1:]):
        linker_mid = d1 + (d2 - d1) // 2
        linker_start = linker_mid - len(linker_params) // 2 - 1
        params[linker_start: linker_start + len(linker_params)] = linker_params
    dna = HelixPose(params)
    return dna


def insert_nucs(dna, dyads, pdb='1KX5'):
    '''
    Fix basepair parameters for nucleosomes in DNA at dyad positions

    Parameters
    ----------
    dna: HelixMC pose
    dyads : dyad positions in DNA (bp)
    Returns
    -------
    dna : HelixMC pose
    nucl : NucleosomeMC pose
    '''
    #    # Get nucleosome
    nucl = nMC.NucPose()
    nucl.from_file(fileio.change_extension(pdb, '3DNA'))

    params = dna.params
    length = len(nucl.dna.params)
    new_dyads = []
    for dyad in dyads:
        start = dyad - nucl.dyad
        end = start + length
        if (start >= 0 and (end < len(params))):
            params[start:end] = nucl.dna.params
            new_dyads.append(dyad)
    dna = HelixPose(params)
    return dna, nucl, np.asarray(new_dyads)


def create_nuc_array(p=None):
    '''
    Create DNA configuration, including nucleosomes

    '''
    if p is None:
        p = Parameters()
        p.add('L_bp', value=500)
        p.add('n_nuc', value=2)
        p.add('NRL', value=197)
    p = p.valuesdict()

    n_bp = p['L_bp']
    n_nucs = p['n_nuc']
    NRL = p['NRL']

    dyads = np.asarray(NRL * (np.arange(0, n_nucs, 1) - (n_nucs - 1) / 2.0))
    dyads = (dyads + n_bp // 2).astype(int)

    random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')
    params = np.tile(random_step.params_avg, (n_bp - 1, 1))

    dna = HelixPose(params)
    dna, nuc, dyads = insert_nucs(dna, dyads)

    return dna, dyads, nuc


def get_fiber_frames(par):
    '''
    Get the frame coordinates that define a stacked, folded fiber

    Returns
    -------
    n_ofs : ndarray of (N, 4, 3)
        The coordinates of the nucleosome frames
    d_ofs : ndarray of (N, 4, 3)
        The coordinates of the dyad frames
    '''

    r = par['diameter_A'].value / 2.0 - 65.0
    nld = par['nld_A'].value
    rise = par['rise_A'].value
    coords = []
    n_ofs = []
    phi = np.sqrt(rise ** 2 - nld ** 2) / r
    if par['chirality'].value < 0:
        phi *= -1
    for i in np.arange(par['n_nuc'].value + 1):
        coords.append(np.asarray([r * np.cos(i * phi), r * np.sin(i * phi), i * nld]))
    for o1, o2 in zip(coords, coords[1:]):
        f0 = par['face'].value * o1 * [-1, -1, 0]
        f0 /= np.linalg.norm(f0)
        f1 = o2 - o1
        f1 /= np.linalg.norm(f1)
        f1 = np.cross(f0, f1)
        f2 = np.cross(f0, f1)
        n_ofs.append(nMC.join_o_f(o1, np.transpose([f0, f1, f2])))
    return n_ofs


def create_folded_fiber(par, nucl):
    dna, dyads, _ = create_nuc_array(p=par)
    params = np.zeros_like(dna.params)
    w = np.zeros(len(dna.coords))
    n_ofs = get_fiber_frames(par)
    origin_of = np.concatenate(([np.zeros(3)], np.eye(3)))
    unwrapped = int(par['Unwrapped_bp'])
    length = len(nucl.dna.params) - 2 * unwrapped
    for n_of, dyad in zip(n_ofs, dyads):
        tf = nMC.get_transformation(nucl.of, target=n_of)
        start_bp = dyad - nucl.dyad + unwrapped
        w[start_bp:start_bp + length] = np.ones(length)
        params[start_bp:start_bp + length] = nucl.dna.params[unwrapped: unwrapped + length]
        params[start_bp - 1] = \
            nMC.ofs2params(origin_of, nMC.apply_transformation(nMC.get_of(nucl.dna, unwrapped), tf))
        params[start_bp + length] = \
            nMC.ofs2params(nMC.apply_transformation(nMC.get_of(nucl.dna, len(nucl.dna.params) - unwrapped), tf),
                           origin_of)
    fiber_dna = HelixPose(params)

    w = np.transpose(np.asarray([w, w, w]))
    return fiber_dna, dyads, w


def residual(fiber_par, ref_of, cast_coords, w):
    coords, _, _ = create_unfolded_fiber(fiber_par, ref_of)
    res = (coords - cast_coords) * w
    sys.stdout.write('.')
    return res


def get_stack_pars(dna, dyads):
    ofs = []
    for dyad in dyads:
        ofs.append(nMC.get_of(dna, dyad))
    params = []
    for of1, of2 in zip(ofs, ofs[1:]):
        params.append(nMC.ofs2params(of1, of2))
    return np.asarray(params)


def create_unfolded_fiber(pars, ref_of):
    dna, dyads, nucl = create_nuc_array(p=pars)
    linker_params = create_curved_linker(p=pars)
    dna = replace_linker(dna, dyads, linker_params)

    if ref_of is None:
        ref_of = nMC.get_of(dna, 0)

    i = pars['ref_frame'].value

    tf = nMC.get_transformation(nMC.get_of(dna, i), target=ref_of)
    coords = nMC.apply_transformation(dna.coords, tf)
    return coords, dna, tf


def scan_fiber_param(fiber_pars, iter_param, iter_vals):
    setnr = 0
    image_files = []
    report_file = fileio.get_filename(incr=True)
    fileio.report_progress(len(iter_vals), title='scan_fiber_param', init=True)
    for i in iter_vals:
        fileio.report_progress(setnr, title='scan_fiber_param\n')

        filename = fileio.get_filename(incr=True, sub=True)
        fiber_pars[iter_param].value = i

        dna, dyads, nucl = create_nuc_array(p=fiber_pars)
        fiber_dna, dyads, w = create_folded_fiber(fiber_pars, nucl)
        fiber_pars['ref_frame'].value = dyads[0]
        ref_of = nMC.get_of(fiber_dna, fiber_pars['ref_frame'].value)

        try:
            out = minimize(residual, fiber_pars, args=(ref_of, fiber_dna.coords, w), method='nelder')
            report_fit(out)
            fiber_pars = out.params
            coords, dna, tf = create_unfolded_fiber(fiber_pars, ref_of)
            image_files.append(
                fileio.pov_dna(filename, [dna.coords], range_A=[400, 400], offset_A=[0, 0, 150], show=True))
            fileio.write_xlsx(filename, str(setnr), fiber_pars, report_file=report_file)
            dna.write2disk(fileio.get_filename(sub=True, ext='npz'))
        except:
            print('Fit did not converge')
        setnr += 1
    fileio.create_mp4_images(image_files, filename=report_file)
    return


def main():
    fiber_par = Parameters()
    fiber_par.add('L_bp', value=400, vary=False)
    fiber_par.add('n_nuc', value=2, vary=False)
    fiber_par.add('nld_A', value=20, vary=False)
    fiber_par.add('NRL', value=197, vary=False)

    fiber_par.add('Unwrapped_bp', value=32, min=0, max=60, vary=False)

    fiber_par.add('Pitch_bp', value=10.380, min=9, max=11, vary=True)
    fiber_par.add('G_kT', value=71.891, min=0, max=100, vary=True)
    fiber_par.add('Phase_rad', value=3.757, min=0, max=2 * np.pi, vary=True)
    fiber_par.add('Period', value=0.787, min=0, max=2, vary=True)
    fiber_par.add('RelFreq', value=1.0, min=0.95, max=1.05, vary=True)

    fiber_par.add('diameter_A', value=330, vary=False)
    fiber_par.add('rise_A', value=100, vary=False)
    fiber_par.add('nld_A', value=30, vary=False)
    fiber_par.add('chirality', value=-1, vary=False)
    fiber_par.add('face', value=1, vary=False)

    fiber_par.add('ref_frame', value=0, vary=False)

    # fiber_par = fileio.read_xlsx('E:\\tmp\\TestData\\data05.xlsx', '16', fiber_par)

    iter_vals = np.linspace(20, 100, 11)
    scan_fiber_param(fiber_par, 'nld_A', iter_vals)
    return


if __name__ == "__main__":
    # execute only if run as a script
    main()
