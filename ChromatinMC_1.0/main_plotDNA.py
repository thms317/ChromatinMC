import FiberMC as FMC


l_handle = 0
n_nuc = 2
NRL = 197
# dna_file = 'simulations\\MC Metropolis-Hastings\\20180220_197_100_NLD.npz'
# dna_file = 'simulations\\NLD\\20180216_197_100_NLD.npz'

dna, dyads, nucl = FMC.create_fiber2(l_handle, n_nuc, NRL, dna_file ='simulations\\MC Metropolis-Hastings\\20180220_197_20_NLD.npz', nuc_file = '1KX5.3DNA', CL = False)
FMC.plotdna2(dna, dyads, nucl)