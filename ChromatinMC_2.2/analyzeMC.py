import matplotlib as mpl

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass
from matplotlib import pyplot as plt
from helixmc.pose import HelixPose
import numpy as np

# ChromatinMC modules:
import NucleosomeMC3 as nMC
import POVutils as pov
import FileIO3 as fileio

default_step_file = 'C:\\Python27\\Lib\\site-packages\\helixmc\\data\\DNA_gau.npy'
default_folder = 'E:\\users\\'
kT = 41



def plot_energy(filename):
    params= np.load(default_step_file)
    p0 = params[0]
    cov = params[1:]
    k = np.linalg.inv(cov)

    sets, files, _ = fileio.contents_xlsx(filename)
    dna = HelixPose.from_file(fileio.change_extension(files[0], 'npz'))

    energy_kT = np.empty((len(files), len(dna.params)))
    i = 0
    fileio.report_progress(len(files), title='analyze_step_parameters', init=True)
    for f in files:
        dna = HelixPose.from_file(fileio.change_extension(f, 'npz'))
        fileio.report_progress(i + 1)
        j = 0
        for p in dna.params:
            energy_kT[i, j] = (np.sum(0.5 * (p - p0)* np.dot(k, p - p0) ))
            j += 1
        i += 1

    F_data = fileio.read_param_xlsx(filename, 'F_pN')

    forces = np.linspace(10, 0, 4)
    energy_F = []
    for force in forces:
        selected = (np.abs(F_data-force)<2.5)
        energy_F.append(np.mean(np.abs(energy_kT[selected]), axis=0))
    energy_F = np.asarray(energy_F)
    print (forces)

    i = xrange(len(dna.params))

    plt.close()
    plt.figure(figsize=(12, 3))
    # plt.plot(i, energy_kT)
    plt.plot(i, energy_F.T)
    energy_thermal = np.ones(len(dna.params)) * 3
    plt.plot(i, energy_thermal, color='k', linestyle=':', linewidth=0.8)

    plt.xlim(0, len(dna.params))
    plt.ylim(-1, 26)
    plt.ylabel('G (kT)')
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0], loc='left',
              fontdict={'fontsize': 10})
    plt.xlabel('i (bp)')

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.draw()
    plt.pause(5)
    filename = fileio.change_extension(filename, '_Edna.jpg')
    plt.savefig(filename, dpi=600, format='jpg')

    return


def plot_step_params(filename, dataset, save=False, wait=0, plot_energy=True):
    filename = change_extension(filename, 'xlsx')
    sets, files, _ = contents_xlsx(filename)
    if dataset == -1:
        filename = sorted(files)[-1]
    else:
        filename = files[sets.index(dataset)]
    dna = HelixPose.from_file(change_extension(filename, 'npz'))

    p0 = np.load(default_step_file)[0]
    sigma2 = np.load(default_step_file)[1:]
    k = 1 / sigma2
    energy_kT = []
    for p in dna.params:
        energy_kT.append(0.5 * (p - p0) * np.dot(k, p - p0) / kT)
    energy_kT = np.asarray(energy_kT)
    energy_kT = np.sum(np.abs(energy_kT), axis=1)

    i = xrange(len(dna.params))

    plt.close()
    plt.figure(figsize=(12, 4))
    if plot_energy:
        plt.plot(i, energy_kT)
        plt.ylim(-1, 25)
        plt.ylabel('G (kT)')
    else:
        plt.plot(i, dna.params)
        plt.ylim(-2.5, 4.5)
        plt.ylabel('step parameters')
        plt.legend(['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist'], loc=1)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0], loc='left',
              fontdict={'fontsize': 10})
    plt.xlabel('i (bp)')

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.draw()
    if wait <> 0:
        plt.pause(wait)
    if save:
        filename = change_extension(filename, '_step.jpg')
        plt.savefig(filename, dpi=600, format='jpg')
    return


def main():
    filename = 'E:\\users\\noort\\data\\20180321\\12nucs_003.dat'
    filename = 'E:\\users\\noort\\data\\20180327\\4x167_001.dat'
    plot_energy(filename)
    return


if __name__ == "__main__":
    # execute only if run as a script
    main()
    # create_pov_mp4('E:\\Users\\Noort\\data\\20180315\\2nucs_001', origin_frame=0, reverse=False)
    # filename = get_filename(incr=True, base='2nucs')
    # print(filename)
