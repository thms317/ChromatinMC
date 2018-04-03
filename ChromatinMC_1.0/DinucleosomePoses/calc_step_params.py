# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:09:58 2016

@author: Visscher
"""

import numpy as np
from matplotlib import pyplot as plt
from helixmc.pose import HelixPose

def plot_params(posename):
     params = np.load(posename)['params']
     plt.close()
     bps = np.arange(len(params))
     parnames = ['shift', 'slide', 'rise',  'tilt', 'roll', 'twist']
     fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, sharex=True)
     plt.suptitle(posename.split('.')[0])
     ax1.scatter(bps, params[:,0])
     ax1.set_title(parnames[0])
     ax2.scatter(bps, params[:,1])
     ax2.set_title(parnames[1])
     ax3.scatter(bps, params[:,2])
     ax3.set_title(parnames[2])
     ax4.scatter(bps, params[:,3])
     ax4.set_title(parnames[3])
     ax5.scatter(bps, params[:,4])
     ax5.set_title(parnames[4])
     ax6.scatter(bps, params[:,5])
     ax6.set_title(parnames[5])
     
     
plot_params('4qlcNRL197.npz')
     