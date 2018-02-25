# -*- coding: UTF-8 -*-

#This file is part of HelixMC.
#    Copyright (C) 2013  Fang-Chieh Chou <fcchou@stanford.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import abc
import warnings
import numpy as np
import NucleosomeMC as nuc
from helixmc import kBT


#####Score function#####
class ScoreBase(object):
    '''
    Base class for scoring fucntion, for inheritence only.
    '''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self):
        return

    @abc.abstractmethod
    def __call__(self, pose):
        return


class ScoreExt(ScoreBase):
    '''
    Score function for force-extension along Z-axis.

    Parameters
    ----------
    force : float
        Applied z-direction force to the helix, in pN.

    Attributes
    ----------
    `force` : float
        See Parameters section above.
    '''
    def __init__(self, force):
        self.force = force

    def __call__(self, pose):
        '''
        Score the input pose.

        Parameters
        ----------
        pose : HelixPose
            Input pose for scoring.

        Returns
        -------
        score : float
            Score of the pose.
        '''
        work = -pose.z_terminal * self.force
#        work = -10000 *kBT
#        print "Work: ", work / kBT
        return work


class ScoreTorsionTrap(ScoreBase):
    '''
    Score function for torsional trap.

    Parameters
    ----------
    stiffness : float
        The stiffness of the torsional trap, in pN.Å.
    target_link : float
        Center of the harmonic torsional trap (link), in radians.

    Attributes
    ----------
    `stiffness` : float
    `target_link` : float
        See Parameters section above.
    '''
    def __init__(self, stiffness, target_link):
        self.stiffness = stiffness
        self.target_link = target_link

    def __call__(self, pose):
        '''
        Score the input pose.

        Parameters
        ----------
        pose : HelixPose
            Input pose for scoring.

        Returns
        -------
        score : float
            Score of the pose.
        '''
        if not pose.compute_tw_wr:
            warnings.warn(
                'pose.compute_tw_wr should be set to Ture for repeating'
                'scoring with target link!!!', RuntimeWarning)
        return (
            0.5 * self.stiffness *
            (pose.link_fuller - self.target_link) ** 2)


class ScoreXyTrap(ScoreBase):
    '''
    Score function for xy trap.

    Parameters
    ----------
    stiffness : float
        The stiffness of the xy trap, in pN/Å.

    Attributes
    ----------
    `stiffness` : float
        See Parameters section above.
    '''
    def __init__(self, stiffness):
        self.stiffness = stiffness

    def __call__(self, pose):
        '''
        Score the input pose.

        Parameters
        ----------
        pose : HelixPose
            Input pose for scoring.

        Returns
        -------
        score : float
            Score of the pose.
        '''
        xy = pose.coord_terminal[:2]
        dist_sq = np.sum(xy ** 2)
        return 0.5 * self.stiffness * dist_sq

class ScoreNucWrap(ScoreBase):
    '''
    Score function for nucleosome wrapping.

    Parameters
    ----------
    energy : float
        The wrapping energy per interaction site in kT
    dyads : ndarray of (N)
        bp position of dyads
    fixed_i : ndarray of (14)
        positions of fixed basepairs relative to dyad
    fixed_params : ndarray of (6)
        step parameters of fixed basepairs

    Attributes
    ----------
    `k` : ndarray of 6
        stiffness for translation and rotation in pN/Å  and 10*kT
    '''
    def __init__(self, wr_energy, dyads, fixed_i, fixed_params):
        self.wr_energy = wr_energy
        self.dyads = dyads
        self.fixed_i = np.asarray(fixed_i)
        self.fixed_params = fixed_params
        
        sd_pos = 0.001 # (A)
        k_pos = 41 / sd_pos**2
        sd_rot = 0.005 #  degrees
        sd_rot = sd_rot *2 * np.pi / 360 # radians
        k_rot = 0.41 / sd_rot*2        
       
        self.k = np.asarray([k_pos, k_pos, k_pos, k_rot, k_rot, k_rot])

    def __call__(self, pose):

        '''
        Score the input pose.

        Parameters
        ----------
        pose : HelixPose
            Input pose for scoring.

        Returns
        -------
        score : float
            Score of the pose.
        '''
        g = 0
        for dyad in self.dyads:
            p = nuc.get_nuc_parameters(pose, self.fixed_i, dyad)
            G = (self.k*(p - self.fixed_params)**2)
            Gtot = np.sum(G, axis = 1)
            for i in range(13):
                Gtot[i] = np.min((Gtot[i], 5*kBT))
#            print Gtot[0:6]/kBT
#            print G
            g = np.sum(Gtot)
#        g = 10000* kBT
#            g=0
        return g


class ScoreAgg(ScoreBase):
    '''
    Score function aggregates of multiple score terms.

    Parameters
    ----------
    score_list : list, optional
        List of score terms (subclass of ScoreBase) in this score.

    Attributes
    ----------
    `score_list` : list
        See Parameters section above.
    `is_empty` : bool
        If the score_list is empty.
    '''
    def __init__(self, score_list=[]):
        self.score_list = score_list

    def __call__(self, pose):
        '''
        Score the input pose.

        Parameters
        ----------
        pose : HelixPose
            Input pose for scoring.

        Returns
        -------
        score : float
            Score of the pose.
        '''
        score = 0
        for term in self.score_list:
            score += term(pose)
        return score

    def append(self, score_term):
        '''
        Append new score term.

        Parameters
        ----------
        score_term : subclass of ScoreBase
            Score term to be appended.
        '''
        self.score_list.append(score_term)

    def clear(self):
        '''
        Clear the score_list.
        '''
        self.score_list = []

    @property
    def is_empty(self):
        return (not self.score_list)


class ScoreTweezers(ScoreAgg):
    '''
    Score function for tweezers experiments.

    Parameters
    ----------
    force : float, optional
        Applied z-stretching force, in pN.
    torsional_stiffness : float, optional
        Stiffness of the torsional trap, in pN.Å.
    target_link : float, optional
        Center of the torsional trap (link), in radians.
    xy_stiffness : float, optional
        Stiffness of xy trap, in pN/Å.
    '''
    def __init__(
            self,
            force=0,
            torsional_stiffness=0,
            target_link=None,
            xy_stiffness=0,
            dyads=[],
            wr_energy = 0,
            fixed_i = [],
            fixed_params = []
    ):
        self.score_list = []
        if force != 0:
            self.score_list.append(ScoreExt(force))
        if torsional_stiffness != 0 and target_link is not None:
            self.score_list.append(
                ScoreTorsionTrap(torsional_stiffness, target_link))
        if xy_stiffness != 0:
            self.score_list.append(ScoreXyTrap(xy_stiffness))
##        if len(dyads) > 1:
##            self.score_list.append(ScoreNucStack(st_energy, dyads, nucleosome))
#        if len(dyads) > 0:
##        if 1 > 0:
#            self.score_list.append(ScoreNucWrap(wr_energy, dyads, fixed_i, fixed_params))
