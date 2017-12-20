from rewotes import crystal, qe, kpointconvg
s = crystal.System()
fe = s.add_species(name='Fe', mass=55.845, potential='Fe.blyp-sp-hgh.UPF')
o = s.add_species(name='O', mass=15.999, potential='O.blyp-hgh.UPF')
s.read_poscar('mp-18731.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
