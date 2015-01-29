# MEPRDC
This repository is a container for scripts (for "Almost" http://open-almost.org/) and configuration files (for "phaistos" http://www.phaistos.org/) which demonstrates the use of Residual Dipolar Coupling data in molecular simulations in a manner consistent with the Maximum Entropy principle. There's also a collective variable for PLUMED2 for use with multiple MD frameworks including gromacs, amber and namd. The exact details of the method is described in a manuscript submitted for publication. A preprint should appear here ASAP. 

Please note that this respository is currently being subject to change, however, previous revisions are stored for reference.

## Prerequisites  
A working installation of almost (from sourceforce branch 2.2), it should be compiled with MPI and LAPACK
Check out the latest version from the Sourceforge repository:
```
svn checkout http://sourceforge.net/p/almost/code/HEAD/tree/branches/almost-2.2 almost

cd almost/

make CXX=mpiCC CC=mpicc CXXFLAGS="-O3 -ffast-math -fomit-frame-pointer -DALM_MPI -DWITH_LAPACK -llapack"

```
_OR_

A working installation of phaistos (we only recommend this option for people who are already familiar with phaistos). Check out the latest version from the Sourceforge respository:
```
svn checkout http://svn.code.sf.net/p/phaistos/code/trunk phaistos
```
note that the _inference_ module should be manually enabled via CMAKE, see phaistos manual for further details.
Please ensure that you have the modules enable which are required by the particular configuration file you wish to use - this information is available in the configuration file of interest. 

_OR_

A working installation of a MD framework supported by PLUMED2 and linked to it. only gromacs is tested at this stage. Please see plumed/ folder for further information

## Future plans

- OpenMP/MPI hybrid parrallelization for the PLUMED2 plugin is not tested! but should work.

## Details 
Almost scripts carry the suffix ".z" which referes to the internal scripting language zeta

Input data format for both phaistos and almost is a space or tab separated text file with the following columns:

1 - Alignment condition identifier

2 - Residue index of first atom

3 - Atom name of first atom

4 - Residue index of second atom

5 - Atom name of second atom

6 - value of observed RDC
 

## Current limitations 

The phaistos suite currently only handles proteins and only single chains. The phaistos example does not allow for automated estimation of the degree of alignment, s.



