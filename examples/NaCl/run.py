from rewotes import crystal, qe, kpointconvg
s = crystal.System()
na = s.add_species(name='Na', mass=22.990, potential='Na.blyp-sp-hgh.UPF')
cl = s.add_species(name='Cl', mass=35.450, potential='Cl.blyp-hgh.UPF')
s.read_poscar('mp-22862.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
