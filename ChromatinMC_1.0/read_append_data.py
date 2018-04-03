import numpy as np
import os

def read_append_data(files_location,file_destination):

    files = os.listdir(files_location)

    array = np.array([])

    for n,f in enumerate(files):
        array = np.concatenate((array, np.loadtxt(files_location+f)))
        n = n

    array = np.transpose(np.hsplit(array,n+1))

    np.savetxt(file_destination+"appended_data.dat",array)

    return

files_location = "C:\\Users\\Brouwer\\Desktop\\ChromatinMC - for PyCharm2\\simulations\\NLD\\dats\\"
file_destination = "C:\\Users\\Brouwer\\Desktop\\ChromatinMC - for PyCharm2\\simulations\\NLD\\"

read_append_data(files_location,file_destination)