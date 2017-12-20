from rewotes import crystal, qe, kpointconvg
s = crystal.System()
ga = s.add_species(name='Ga', mass=69.723, potential='Ga.blyp-d-hgh.UPF')
as_ = s.add_species(name='As', mass=74.921, potential='As.blyp-hgh.UPF')
s.read_poscar('mp-10048.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
