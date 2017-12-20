from rewotes import crystal, qe, kpointconvg
s = crystal.System()
si = s.add_species(name='Si', mass=28.0855, potential='Si.blyp-hgh.UPF')
s.read_poscar('mp-10649.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
