# TAPPS
Thermophysical Ab-initio Properties from Phase Shifts

This repository contains a set of Python3 programs and corresponding data files enabling the calculation of thermophysical properties of noble gases from phase shift and bound state data.

The main files are:

* ``thermophysicalPairProperties.py`` to compute the second virial coefficient, its first two temperature derivatives and the second acoustic virial coefficient for the pure and cross phases.
* ``Ne_avg_virials.py`` to compute the second virial coefficient, its first two temperature derivatives and the second acoustic virial coefficient for the normal mixture of neon isotopes.
* ``data`` directory with compressed JSON files containing phase shifts and bound states computed for
  - ${}^3$He, ${}^4$He, ${}^3$He-${}^4$He using the pair potential from [Czachorowski et al.](https://doi.org/10.1103/PhysRevA.102.042810)
  - ${}^{20}$Ne, ${}^{21}$Ne, ${}^{22}Ne and their cross values, using the pair potential from [Hellmann et al.](https://doi.org/10.1063/5.0047999)
  - ${}^{40}$Ar using the pair potential from [Lang et al.](https://doi.org/10.1103/PhysRevA.109.052803)

The ``README.TXT`` file contains detailed instructions on how to use the Python programs and the classes contained therein.
