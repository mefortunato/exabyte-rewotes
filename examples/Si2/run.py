from rewotes.kpointconvg import System, KPointConvg
s = System()
si = s.add_species(name='Si', mass=28.0855, potential='si_pbe_gbrv_1.0.upf')
s.read_poscar('si.poscar')
kconvg = KPointConvg(system=s)
kdim = kconvg.find_convg(tol=1e-5)

print('Total energy converged with tolerace 1e-5 eV at k dimensions of {}'.format(kdim))
