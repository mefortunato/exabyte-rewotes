from rewotes import crystal, qe, kpointconvg
s = crystal.System()
mg = s.add_species(name='Mg', mass=24.305, potential='Mg.blyp-hgh.UPF')
o = s.add_species(name='O', mass=15.999, potential='O.blyp-hgh.UPF')
s.read_poscar('mp-1265.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
