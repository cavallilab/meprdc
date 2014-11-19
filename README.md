# MEPRDC
This repository is a container for scripts (for "Almost" http://open-almost.org/) and configuration files (for "phaistos" http://www.phaistos.org/) which demonstrates the use of Residual Dipolar Coupling data in molecular simulations in a manner consistent with the Maximum Entropy principle. The exact details of the method is described in a manuscript submitted for publication. A preprint should appear here ASAP. 

Please note that this respository is currently being subject to change, however, previous revisions will be stored for reference.

## Prerequisites  
A working installation of almost (from sourceforce branch 2.2)
Check out the latest version from the Sourceforge repository:
```
svn checkout http://sourceforge.net/p/almost/code/HEAD/tree/branches/almost-2.2 almost
```

_OR_

A working installation of phaistos. Check out the latest version from the Sourceforge respository:
```
svn checkout http://svn.code.sf.net/p/phaistos/code/trunk phaistos
```
note that the _inference_ module should be manually enabled via CMAKE, see phaistos manual for further details.
Please ensure that you have the modules enable which are required by the particular configuration file you wish to use - this information is available in the configuration file of interest. 

## Future plans
- Implementation of a collective variable for PLUMED2

## Details 
Almost scripts carry the suffix ".z" which referes to the internal scripting language zeta

 

## Current limitations 

The phaistos suite currently only handles proteins and only single chains. The phaistos example does not allow for automated estimation of the degree of alignment, s.



