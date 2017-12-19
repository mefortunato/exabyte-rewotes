from __future__ import print_function

import os
import sys
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE


class KPointConvgException(Exception):
    pass


class System(object):
    """rewotes.kpointconvg.System

    System class defining a crystal structure for DFT calculation.
    
    Attributes:
        species: atomic species list of dict objects defining {'name': str, 'mass': float, 'potential': str}
        atoms: list of atoms dict objects defining {'species': species, 'coordinates': ndarray(shape=(3,))}
        cell: ndarray(shape=(3,)) of cell vectors
    """
    def __init__(self, species=None, atoms=None, cell=None):
        self.species = species
        self.atoms = atoms
        self.cell = cell
        
        if self.species is None:
            self.species = []
        if self.atoms is None:
            self.atoms = []
        if self.cell is None:
            self.cell = np.zeros(shape=(3,3))
        
    def read_poscar(self, fname):
        """rewotes.kpointconvg.System.read_poscar

        Iterates through s.bonds, s.angles, s.dihedrals, and s.impropers and removes
        those which contain this :class:`~pysimm.system.Particle`.

        Args:
            s: :class:`~pysimm.system.System` object from which bonding objects will be removed

        Returns:
            None
        """
        with open(fname) as f:
            for _ in range(2):
                line = f.next()
            self.cell = np.zeros(shape=(3,3))
            for n in range(3):
                line = f.next()
                self.cell[n] = map(float, line.split())
            line = f.next()
            line = f.next()
            natoms = np.array(map(int, line.split())).sum()
            line = f.next()
            for n in range(natoms):
                line = f.next().split()
                species = None
                for sp in self.species:
                    if sp['name'] == line[-1]:
                        species = sp
                        break
                if species is None:
                    print('cannot find species {} in system'.format(line[-1]))
                    return None
                self.add_atom(species=species, coordinates=map(float, line[:3]))
                
        
    def add_species(self, name, mass, potential):
        if not os.path.isfile(os.path.join('pseudo', potential)):
            raise KPointConvgException('Cannot find pseudopotential file {} in pseudo dir'.format(potential))
        species = {
            'name': name,
            'mass': mass,
            'potential': potential
        }
        self.species.append(species)
        return species
        
    def add_atom(self, species, coordinates):
        atom = {
            'species': species,
            'x': coordinates[0],
            'y': coordinates[1],
            'z': coordinates[2]
        }
        self.atoms.append(atom)
        return atom

    def write_input(self):
        s = 'ATOMIC_SPECIES\n'
        for sp in self.species:
            s += '{name} {mass} {potential}\n'.format(**sp)
        s += 'ATOMIC_POSITIONS crystal\n'
        for a in self.atoms:
            s += '{species_name} {x} {y} {z}\n'.format(species_name=a['species']['name'], **a)
        s += 'CELL_PARAMETERS angstrom\n'
        for v in self.cell:
            s += '{} {} {}\n'.format(*v)
        return s
        

class KPointConvg(object):
    def __init__(self, system=None, kmin=1, kmax=10, PW_EXEC='pw.x', kinetic_cutoff=40):
        self.system = system
        self.kmin = kmin
        self.kmax = kmax
        self.kinetic_cutoff = kinetic_cutoff
        self.PW_EXEC = PW_EXEC
        self.data = pd.DataFrame(columns=['k1', 'k2', 'k3'])
        
    @property
    def _input(self):
        s = '''&CONTROL
    calculation = 'scf'
    title = ''
    verbosity = 'low'
    restart_mode = 'from_scratch'
    wf_collect = .true.
    tstress = .true.
    tprnfor = .true.
    outdir = 'outdir'
    wfcdir = 'outdir'
    prefix = '__prefix__'
    pseudo_dir = 'pseudo'
/
'''
        s += '''&SYSTEM
    ibrav = 0
    nat = {}
    ntyp = {}
    ecutwfc = {}
    ecutrho = 200
    occupations = \'smearing\'
    degauss = 0.005
/
'''.format(len(self.system.atoms), len(self.system.species), self.kinetic_cutoff)

        s += '''&ELECTRONS
    diagonalization = 'david'
    diago_david_ndim = 4
    diago_full_acc = .true.
    mixing_beta = 0.3
    startingwfc = 'atomic+random'
/
'''
        s += self.system.write_input()
        s += '''K_POINTS automatic
{} {} {} 0 0 0
'''
        return s
        
    def run(self, k, input_file=None, debug=False, nproc=1, mpi_prefix='mpirun', save_input='tmp.in', save_output='tmp.out'):
        if type(k) == int:
            k1 = k2 = k3 = k
        elif type(k) == list or (type(k) == np.ndarray and k.shape == (3,)):
            k1, k2, k3 = k
        else:
            print('cannot understand kpoints')
            return None
        if input_file is None:
            inp = self._input.format(k1, k2, k3)
            if debug:
                print(inp)
            with open(save_input, 'w') as f:
                f.write(inp)
            input_file = save_input
        if nproc > 1:
            cmd = [mpi_prefix, '-np', str(nproc), self.PW_EXEC, '-inp', input_file]
        else:
            cmd = [self.PW_EXEC, '-inp', input_file]
        print('running with k=[{} {} {}]...'.format(k1, k2, k3))
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        stdo, stde = p.communicate()
        
        with open(save_output, 'w') as f:
            f.write(stdo)
            f.write(stde)
        
        if debug:
            print(stdo)
            print(stde)
        
        lines = stdo.split('\n')
        E_Ry = None
        for line in lines:
            if line.startswith('!    total energy'):
                E_Ry = float(line.split('=')[1].split('Ry')[0].strip())
        if E_Ry is not None:
            self.data = self.data.append({'k1': k1, 'k2': k2, 'k3': k3, 'E_Ry': E_Ry, 'E_eV': E_Ry/13.605698066}, ignore_index=True)
            
    def find_convg(self, criteria='E_eV', tol=0.00001, debug=False, nproc=1, mpi_prefix='mpiexec'):
        for i in range(self.kmin, self.kmax+1):
            self.run(i, debug=debug, nproc=nproc, mpi_prefix=mpi_prefix, save_input='{}.{}.{}.in'.format(i, i, i), save_output='{}.{}.{}.out'.format(i, i, i))
            if len(self.data) > 1:
                self.data['delta_{}'.format(criteria)] = self.data.diff()[criteria].abs()
                if np.any(self.data.loc[self.data['delta_{}'.format(criteria)] < tol]):
                    return list(self.data.loc[self.data['delta_{}'.format(criteria)] < tol, ['k1', 'k2', 'k3']].values[0])
        print('Did not converge')
        return None