K-Point Convergence Tracker
===========================

This is a solution to one of Exabyte's REal WOrld TEstS (rewotes) described [here](https://github.com/Exabyte-io/rewotes/blob/master/Convergence-Tracker.md).
In summary, the code is supposed to find a convergence criteria as a function of the k-point parameter for a DFT calculation.  

This code has been designed to be used with the open-source code for electronic structure calculations, [Quantum ESPRESSO](http://www.quantum-espresso.org/).
The execution of a calculation using Quantum ESPRESSO is handled by the `qe` module. Extension for this kpoint tracker to other calculation engines would require development of a similar module for other software.

Dependencies
============

This code depends on the proper installation of the [Quantum ESPRESSO](http://www.quantum-espresso.org/) software package. See their [installation directions](http://www.quantum-espresso.org/download/) for further information on installation.

The code has been tested for both Python2.7.14 and Python3.5.2.

The following Python packages are dependencies and must be installed prior to use of this package:
* [NumPy](http://www.numpy.org/)
* [Pandas](http://pandas.pydata.org/)

Usage
=====

There are three modules that should be imported: `crystal`, `qe`, and `kpointconvg`:

```python
from rewotes import crystal, qe, kpointconvg
```

Begin by defining a crystal structure using the `crystal.System` class and adding the atomic species by defining a name, mass, and pseudopotential file:
```python
system = crystal.System()
si = system.add_species(name='Si', mass=28.0855, potential='si_pbe_gbrv_1.0.upf')
```
**NOTE:** The pseudopotential file should be located in a directory named 'pseudo' in the working directory where you will run your calculation.
An exeption will be raised if pseudo/{species['potential']} does not exist.

Atoms can be added to the system either manually:
```python
a1 = system.add_atom(species=si, coordinates=[0., 0., 0.])
a2 = system.add_atom(species=si, coordinates=[0.25, 0.25, 0.25])
```

or by reading a poscar file (obtained from Exabyte Materials builder):
```python
system.read_poscar('si.poscar')
```

Then build a clalculation using the `qe.Calculation` class, passing a `crystal.system` object and calculation settings:
```python
calc = qe.Calculation(system, kinetic_cutoff=40)
```

The `kpointconvg.KPointConvg` class should then be instantiated by passing the `Calculation` object, and optionally kmin and kmax to search for convergence (defaults: kmin=1, kmax=10):
```python
kconvg = kpointconvg.KPointConvg(calc)
```

Call the find_convg() method to search for convergence optionally passing a tolerance (default: 1e-5 eV):
```python
kconvg.find_convg(tol=1e-5)
```

The find_cong() method will return either the kpoint dimensions that satisfy the convergence criteria, or will print 'Did not converge' and return `None`

Advanced Usage
==============

Currently only the total energy and pressure is extracted from the Quantum ESPRESSO output. To extend this package to find the convergence with respect to other properties (such as forces, etc.), functionality to extract those properties from the output must be implemented in the `qe` module.
For this implementation, this data should be included in the Pandas DataFrame (stored as KPointConvg.data), in which case the user can specify this new convergence criteria when calling the find_convg() method with the keyword 'criteria' (which by default is 'E_eV, the name of the column in KPointConvg.data).

To find convergence with respect to pressure within 1kbar:

```python
kconvg.find_convg(critera='pressure', tol=1.)
```