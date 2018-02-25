# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 16:19:13 2016

@author: Visscher
"""
'''
simulation dates

3kT 1000 runs
(14,53)
20160414

4kT 1000 runs
(1,12)
20160419

20kT 50 runs
(15,30)
20160420

3kT 50 runs
(1,16)
20160421

10kT 50 runs
(17,32)
20160421

3kT 1000 runs 
(37,38)

'''

#Read data files from the directory and put the result into a single file
#Creat

import numpy as np
from matplotlib import pyplot as plt
import os


kBT = 4.114 #pn*nm

def WLC(L):
    z = np.linspace(0, .95*L)
    A = 50.0 #nm
    f = kBT/A*(.25*(1.0/(1-z/L)**2-1.0)+z/L)
    return z, f

data_dir = "C:\\Users\\Visscher\data\\"

date_dir = "20160421\\"


start = 15
end = 30
data = []
for i in range(start, end+1):
    file_name = 'data_{:0>4}.dat'.format(str(i))
    fpath = os.path.join(data_dir+date_dir, file_name)
    data.append(np.loadtxt(fpath, skiprows = 1))
    
data = np.transpose(np.array(data))
name = 'compare_points'
np.savetxt(name + '.dat', data, header = 'Force mean(Z) std(Z)')







#plt.plot(data[1],data[0], 'ro', label = '20kT')
#plt.plot(data1[1],data1[0], 'bo', label = "3kT")
#plt.plot(data2[1],data2[0], 'go', label = "10kT")
#plt.legend(loc = 0)
#plt.xlabel('Z (nm)')
#plt.ylabel('Force(nm)')
#plt.title(name)
#plt.show()
