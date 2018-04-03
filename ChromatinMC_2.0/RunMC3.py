# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 10:07:46 2018

@author: noort
"""
from __future__ import print_function
import matplotlib as mpl

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass

import numpy as np
from lmfit import Parameters
from helixmc import util
from helixmc.random_step import RandomStepSimple, RandomStepAgg, symmetrize_WC
from helixmc.score import ScoreTweezers
from helixmc.pose import HelixPose
# ChromatinMC modules:
import FiberMC5 as fMC
import NucleosomeMC3 as nMC
import FileIO3 as fileio

kT = 41.0
np.set_printoptions(formatter={'float': '{: 0.1f}, '.format})


def get_unwrap_energy(wrap_params, fixed_wrap_params, g_wrap_kT):
    sigma_trans = 1  # 1 Angstrom, similar to basepair steps
    sigma_rot = 5 * np.pi / 180  # 5 degree, similar to basepair steps
    k1 = kT / np.asarray([sigma_trans, sigma_trans, sigma_trans, sigma_rot, sigma_rot, sigma_rot]) ** 2
    G_unwrap = 0.5 * k1 * (wrap_params - fixed_wrap_params) ** 2
    G_unwrap = np.clip(G_unwrap, 0, g_wrap_kT * kT / 6)
    # print np.sum(G_unwrap, axis=1)/kT
    return np.sum(G_unwrap), G_unwrap


def score_wrapping(moving_bp, dna, dyads, nucl, fixed_wrap_params, g_wrap_kT):
    closest_dyad = dyads[np.abs(dyads - moving_bp).argmin()]
    start_bp = closest_dyad - nucl.dyad
    end_bp = start_bp + len(nucl.dna.params)
    if start_bp <= moving_bp < end_bp:
        wrap_params = nMC.get_wrap_params(dna, closest_dyad, nucl.fixed)
        G, _ = get_unwrap_energy(wrap_params, fixed_wrap_params, g_wrap_kT)
        return G
    else:
        return 0


def score_stacking(moving_bp, dna, dyads, fixed_stack_params, g_stack_kT, fiber_start=1):
    start_dyad = np.argmax(dyads > moving_bp) - 1
    if 0 <= start_dyad < len(dyads) - fiber_start:
        stack_params = fMC.get_stack_pars(dna, [dyads[start_dyad], dyads[start_dyad + fiber_start]])
        sigma = np.asarray([1, 1, 1, 0.1, 0.1, 0.1])
        sigma /= 1
        k = kT / sigma ** 2
        g = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g = np.clip(g, 0, g_stack_kT * kT)
        # second potential
        sigma = np.asarray([10, 10, 10, 1e3, 1e3, 1e3])
        k = kT / sigma ** 2
        g2 = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g2 = np.clip(g2, 0, g_stack_kT * kT)
        g += g2
    else:
        g = 0
    return g


def score_surface(dna):
    surface = dna.coords[:, 2] < 0
    if np.sum(surface) > 0:
        return 1e7
    else:
        return 0


def get_new_step_params(moving_bp, prev_bp, dna, dyads, nucl, random_step):
    closest_dyad = dyads[np.abs(dyads - moving_bp).argmin()]
    nucl_bp = moving_bp - closest_dyad + nucl.dyad
    nucl_prev_bp = prev_bp - closest_dyad + nucl.dyad
    if 1 <= nucl_bp < len(nucl.dna.params) - 1:
        if np.array_equal(dna.params[prev_bp], nucl.dna.params[nucl_prev_bp]):
            new_step_params = nucl.dna.params[nucl_bp]
        else:
            new_step_params = random_step()[0]
    else:
        new_step_params = random_step()[0]
    return new_step_params


def MC_move(dna, moving_bp, previous_bp, scorefxn, fixed_wrap_params, fixed_stack_params, dyads, nucl,
            random_step, g_wrap_kT, g_stack_kT, fiber_start):
    old_step_params = dna.params[moving_bp]
    new_step_params = get_new_step_params(moving_bp, previous_bp, dna, dyads, nucl, random_step)
    if np.array_equal(old_step_params, new_step_params):
        return False
    else:
        old_score = scorefxn(dna) \
                    + score_wrapping(moving_bp, dna, dyads, nucl, fixed_wrap_params, g_wrap_kT) \
                    + score_stacking(moving_bp, dna, dyads, fixed_stack_params, g_stack_kT, fiber_start) \
                    + score_surface(dna)
        dna.update(moving_bp, new_step_params)
        new_score = scorefxn(dna) \
                    + score_wrapping(moving_bp, dna, dyads, nucl, fixed_wrap_params, g_wrap_kT) \
                    + score_stacking(moving_bp, dna, dyads, fixed_stack_params, g_stack_kT, fiber_start) \
                    + score_surface(dna)
        if util.MC_acpt_rej(old_score, new_score):
            return True
        else:
            dna.update(moving_bp, old_step_params)
            return False


def main():
    # Define fiber params
    pars = Parameters()
    pars.add('F_pN', value=0)
    pars.add('z_nm', value=0)
    pars.add('L_bp', value=800)
    pars.add('n_nuc', value=4)
    pars.add('NRL', value=167)
    pars.add('P_nm', value=50)
    pars.add('dyad0_bp', value=0)
    pars.add('diameter_A', value=330)
    pars.add('rise_A', value=100)
    pars.add('nld_A', value=25)
    pars.add('chirality', value=-1)
    pars.add('Unwrapped_bp', value=30)
    pars.add('face', value=1)
    pars.add('g_wrap_kT', value=3)
    pars.add('g_stack_kT', value=25)
    pars.add('fiber_start', value=2)

    g_wrap_kT = pars['g_wrap_kT'].value
    g_stack_kT = pars['g_stack_kT'].value

    # Setup files and forces
    filename = fileio.get_filename(incr=True, root='3nucs')
    print(filename)
    n_step = 250
    n_substeps = 10
    fmax_pN = 10
    fmin_pN = 0.1
    # forces = np.linspace(fmin_pN, fmax_pN, n_step / 2)
    forces = np.logspace(np.log10(fmin_pN), np.log10(fmax_pN), n_step / 2)
    forces = np.append(forces, forces[::-1])

    get_from_file = False
    if get_from_file:
        # Get from file
        file_in = 'E:\\Users\\noort\\data\\20180321\\data_003.xlsx'
        dataset = 0
        pars, datafile = fileio.read_xlsx(file_in, dataset, pars=pars)
        dna, dyads, nucl = fMC.create_nuc_array(p=pars)
        dna = HelixPose.from_file(fileio.change_extension(datafile, 'npz'))
    else:
        # Initialize fiber pose
        dna, dyads, nucl = fMC.create_nuc_array(p=pars)

    pars['dyad0_bp'].value = dyads[0]
    fileio.plot_dna(dna, title='Initial conformation\n', range_nm=100, save=True)
    pars.pretty_print(columns=['value'])

    # Get stack and wrap parameters
    fixed_wrap_params = nMC.get_wrap_params(nucl.dna, nucl.dyad, nucl.fixed)
    fiber_dna, dyads, w = fMC.create_folded_fiber(pars, nucl)
    fixed_stack_params = fMC.get_stack_pars(fiber_dna, dyads)[0]
    fiber_start = pars['fiber_start'].value
    # Initialize random steps
    random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')

    basepairs = np.asarray(xrange(pars['L_bp'] - 1))
    accept = 0
    all_coord = np.empty((n_step, 3))

    current_step = 0
    fileio.report_progress(n_step, title='RunMC3', init=True)
    for force in forces:
        fileio.report_progress(current_step + 1, title='RunMC3')
        scorefxn = ScoreTweezers(force)
        previous_bp = 0
        for sub_step in xrange(n_substeps):
            for bp in basepairs:
                accept += MC_move(dna, bp, previous_bp, scorefxn, fixed_wrap_params, fixed_stack_params,
                                  dyads, nucl, random_step, g_wrap_kT, g_stack_kT, fiber_start)
                previous_bp = bp
            basepairs = basepairs[::-1]

        fileio.plot_dna(dna, update=True, title='F = {:.1f} pN\n'.format(force), save=True)

        pars['F_pN'].value = force
        pars['z_nm'].value = dna.coord_terminal[2]
        fileio.write_xlsx(fileio.get_filename(sub=True), str(current_step), pars, report_file=filename)
        dna.write2disk(fileio.get_filename(sub=True, ext='npz'))
        all_coord[current_step] = dna.coord_terminal
        current_step += 1

    z = all_coord[:, 2] / 10
    wlc = 1 - 0.5 * np.sqrt(0.1 * kT / (forces * pars['P_nm']))
    grid = []
    for i in xrange(1, pars['n_nuc'] + 1):
        grid.append(wlc * (pars['L_bp'] - 80 * i) / 3)
        grid.append(wlc * (pars['L_bp'] - (pars['n_nuc'] * 80 + i * (147 - 80))) / 3)
    wlc *= pars['L_bp'] / 3
    selected = np.diff(np.append([-1], forces)) > 0

    fileio.save_plot((forces, z, wlc, selected), filename=filename,
                     ax_labels=['z (nm)', 'F (pN)'], grid=grid, transpose=True, xrange=[0, 1.1 * pars['L_bp'] / 3])

    if n_step > 0:
        print('Accept rate = %.1f %%' % (100 * float(accept) / (n_step * n_substeps * (pars['L_bp'] - 1))))
    try:
        fileio.create_mp4_pov(fileio.get_filename(sub=True, folder=True), origin_frame=0, reverse=False)
    except Exception, e:
        print(Exception, e)


if __name__ == '__main__':
    main()
