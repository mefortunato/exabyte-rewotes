from rewotes import crystal, qe, kpointconvg
s = crystal.System()
ti = s.add_species(name='Ti', mass=47.867, potential='Ti.blyp-sp-hgh.UPF')
o = s.add_species(name='O', mass=15.999, potential='O.blyp-hgh.UPF')
s.read_poscar('mvc-12939.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
