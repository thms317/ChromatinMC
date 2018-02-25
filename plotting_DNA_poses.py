import numpy as np
import FiberMC as FMC
import matplotlib.pyplot as plt

dna_file = 'simulations\\20180215_197_0_unwrap.npz'

DNA, dyads, nucl = FMC.create_fiber2(0, 2, 197, dna_file=dna_file, nuc_file='1KX5.3DNA', CL=False)
# FMC.plotdna2(DNA, dyads, nucl, plt_range = 500)

unwrap = 0

start = dyads[0] + (147 - nucl.d_index)
end = dyads[1] - nucl.d_index


def GetCut(DNA, nucl, n_nuc, cut, dyads, diameter=330.0, rise=100.0, nld=100, ribs=1):
    fiber_params = FMC.set_fiber_params(n_nuc, diameter=diameter, rise=rise, nld=nld, ribs=ribs)
    # fiber_params[0,5] = 1.6
    fiber_dyad_frames = FMC.get_fiber_dyad_frames(fiber_params, nucl.d_index, nucl.nuc_type)
    FMC.cast_nucs_fiber(DNA, dyads, fiber_dyad_frames, cut=cut)


def MC(DNA, start=179, end=239, unwrap=-20, diameter=330.0, rise=100.0, nld=100, ribs=1):
    kBT = 4.1 * 10  # We work in Angstrom
    for i in range(start + unwrap, end - unwrap):
        Estart = FMC.Boltzmann_Energy(DNA.params[start + unwrap]) + FMC.Boltzmann_Energy(DNA.params[i])
        pars = DNA.params
        pars[i] += (0.5 - np.random.rand()) * np.array([0.001, 0.001, 0.01, 0.001, 0.001, 0.01])  # Bert's steps
        # pars[i] += (0.5 - np.random.rand()) * np.array([0.01, 0.01, 0.1, 0.01, 0.01, 0.1])
        newDNA = DNA.copy()
        newDNA.set_params(pars)

        GetCut(newDNA, nucl, 2, start - dyads[0] + unwrap, dyads, diameter=diameter, rise=rise, nld=nld, ribs=ribs)
        Enew = FMC.Boltzmann_Energy(newDNA.params[start + unwrap]) + FMC.Boltzmann_Energy(newDNA.params[i])

        if Enew - Estart <= 0:
            R = 1
            print Estart - Enew
        elif Enew - Estart >= 0:
            # calculate ratio
            P_t = np.exp(-Enew / kBT)
            P_n = np.exp(-Estart / kBT)
            R = P_t / P_n
            print Estart - Enew

        pars[i] = pars[i] * R
        newDNA = DNA.copy()
        newDNA.set_params(pars)
        DNA.set_params(newDNA.params)


#
#     pars = DNA.params
#     return pars
#
# pars = MC(DNA, start=start, end=end, unwrap=unwrap, diameter=330.0, rise=100.0, nld=100, ribs=1)
#
# pars2 = [list(item) for item in pars]
# pars2 = np.array(pars2)

file_name = 'test'

E_link = []
E = []
E_bp = []
number_of_iterations = 10
j = 0
while j < number_of_iterations:
    MC(DNA, start, end, unwrap)
    # print str(j + 1) + ' of ' + str(number_of_iterations) + ' (' + str(nld) + ' nucleosome line density)'
    # nice to include the energy calculated from the helixmc score-function
    Etemp = FMC.Boltzmann_Energy(DNA.params[start + unwrap])
    print Etemp
    E.append(Etemp)

    E_l = []
    for i in range(start + unwrap, end - unwrap):
        E_l.append(FMC.Boltzmann_Energy(DNA.params[i]))

    E_link.append(np.sum(E_l))
    E_bp = np.concatenate((E_bp, np.array(E_l)))

    j += 1

nld = 0

np.savetxt(str(file_name) + '_' + str(nld) + '_NLD_E.dat ', E)
np.savetxt(str(file_name) + '_' + str(nld) + '_NLD_E_linkerDNA.dat ', E_link)
DNA.write2disk(str(file_name) + '_' + str(nld) + '_NLD')

E_bp = np.transpose(np.hsplit(E_bp, j))
np.savetxt(str(file_name) + '_basepair_energies_NLD' + str(nld) + '.dat', E_bp)

# E_res = np.concatenate((E_res, E))
# E_link_series = np.concatenate((E_link_series, E_link))
