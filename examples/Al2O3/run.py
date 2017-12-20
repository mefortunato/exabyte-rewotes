from rewotes import crystal, qe, kpointconvg
s = crystal.System()
al = s.add_species(name='Al', mass=26.982, potential='')
o = s.add_species(name='O', mass=15.999, potential='Al.blyp-hgh.UPF')
s.read_poscar('mp-642363.poscar')
calc = qe.Calculation(s)
kconvg = kpointconvg.KPointConvg(calc)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
