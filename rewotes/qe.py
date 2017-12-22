from __future__ import print_function

import re
from subprocess import Popen, PIPE


def write_system(system):
    """rewotes.qe.write_system

    Formats string representation of system formatted as input for Quantum ESPRESSO

    Args:
        crystal.System object

    Returns:
        Input string
    """
    s = 'ATOMIC_SPECIES\n'
    for sp in system.species:
        s += '{name} {mass} {potential}\n'.format(**sp)
    s += 'ATOMIC_POSITIONS crystal\n'
    for a in system.atoms:
        s += '{species_name} {x} {y} {z}\n'.format(species_name=a['species']['name'], **a)
    s += 'CELL_PARAMETERS angstrom\n'
    for v in system.cell:
        s += '{} {} {}\n'.format(*v)
    return s


class Calculation(object):
    def __init__(self, system, kinetic_cutoff=40, kpoints=[1, 1, 1], exec_='pw.x', mpi_prefix='mpirun'):
        self.system = system
        self.kinetic_cutoff = kinetic_cutoff
        self.kpoints = kpoints
        self.exec_ = exec_
        self.mpi_prefix = mpi_prefix
        
    @property
    def input(self):
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
        s += write_system(self.system)
        s += '''K_POINTS automatic
{} {} {} 0 0 0
'''.format(*self.kpoints)
        return s
        
    def run(self, input_file=None, debug=False, nproc=1, save_input='tmp.in', save_output='tmp.out'):
        if input_file is None:
            if debug:
                print(self.input)
            with open(save_input, 'w') as f:
                f.write(self.input)
            input_file = save_input
        if nproc > 1:
            cmd = [self.mpi_prefix, '-np', str(nproc), self.exec_, '-inp', input_file]
        else:
            cmd = [self.exec_, '-inp', input_file]
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        stdo, stde = p.communicate()
        
        with open(save_output, 'w') as f:
            f.write(stdo)
            f.write(stde)
        
        if debug:
            print(stdo)
            print(stde)
            
        # TODO
        # catch exception if convergence not achieved
            
        m = re.search(r'!.*total energy.*=(.+?)Ry', stdo)
        try:
            E_Ry = float(m.group(1))
        except:
            E_Ry = None
        
        m = re.search('P=(.+?)\n', stdo)
        try:
            press = float(m.group(1))
        except:
            press = None
            
        return {'E_Ry': E_Ry, 'pressure': press, 'E_eV': E_Ry/13.605698066}