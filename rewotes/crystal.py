from __future__ import print_function

import os
import numpy as np


class CrystalException(Exception):
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

        Read atom species names, atom coorindates and cell dimensions from poscar file format.

        Args:
            fname: path to poscar file

        Returns:
            None
        """
        with open(fname) as f:
            for _ in range(2):
                line = next(f)
            self.cell = np.zeros(shape=(3,3))
            for n in range(3):
                line = next(f)
                self.cell[n] = list(map(float, line.split()))
            line = next(f)
            line = next(f)
            natoms = np.array(list(map(int, line.split()))).sum()
            line = next(f)
            for n in range(natoms):
                line = next(f).split()
                species = None
                for sp in self.species:
                    if sp['name'] == line[-1]:
                        species = sp
                        break
                if species is None:
                    raise CrystalException('Cannot find species {} in system'.format(line[-1]))
                self.add_atom(species=species, coordinates=list(map(float, line[:3])))
                
        
    def add_species(self, name, mass, potential):
        """rewotes.kpointconvg.System.add_species

        Adds atomic species to System.species.

        Args:
            name: species name (str)
            mass: mass of atomic species (float)
            potential: name of pseudopotential file located in directory names pseudo (str)

        Returns:
            Dictionary representing species for reference later
        """
        if not os.path.isfile(os.path.join('pseudo', potential)):
            raise CrystalException('Cannot find pseudopotential file {} in pseudo dir'.format(potential))
        species = {
            'name': name,
            'mass': mass,
            'potential': potential
        }
        self.species.append(species)
        return species
        
    def add_atom(self, species, coordinates):
        """rewotes.kpointconvg.System.add_species

        Adds atom to System.atoms.

        Args:
            species: dictionary representing species (should have been added to system already and contain 'name')
            coordinates: array of length 3 with fractional coordinates of atom

        Returns:
            Dictionary representing atom for reference later
        """
        atom = {
            'species': species,
            'x': coordinates[0],
            'y': coordinates[1],
            'z': coordinates[2]
        }
        self.atoms.append(atom)
        return atom