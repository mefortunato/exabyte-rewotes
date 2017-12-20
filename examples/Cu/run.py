from rewotes import crystal, qe, kpointconvg
s = crystal.System()
fe = s.add_species(name='Cu', mass=63.546, potential='Cu.blyp-d-hgh.UPF')
s.read_poscar('mp-989782.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
