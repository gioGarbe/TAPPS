# TAPPS
**T**hermophysical **A**b-initio **P**roperties from **P**hase **S**hifts

This repository contains a set of Python3 programs and corresponding data files enabling the calculation of thermophysical properties of noble gases from phase shift and bound state data.

## Main files

The main files are:

* ``thermophysicalPairProperties.py`` to compute the second virial coefficient, its first two temperature derivatives and the second acoustic virial coefficient for the pure and cross phases.
* ``Ne_avg_virials.py`` to compute the second virial coefficient, its first two temperature derivatives and the second acoustic virial coefficient for the normal mixture of neon isotopes.
* ``data/`` directory with compressed JSON files containing phase shifts and bound states computed for
  - <sup>3</sup>He, <sup>4</sup>He, <sup>3</sup>He–<sup>4</sup>He using the pair potential from [Czachorowski et al.](https://doi.org/10.1103/PhysRevA.102.042810)
  - <sup>20</sup>Ne, <sup>21</sup>Ne, <sup>22</sup>Ne and their cross values, using the pair potential from [Hellmann et al.](https://doi.org/10.1063/5.0047999)
  - <sup>40</sup>Ar using the pair potential from [Lang et al.](https://doi.org/10.1103/PhysRevA.109.052803)

## Documentation

The ``README.TXT`` file contains detailed instructions on how to use the Python programs and the classes contained therein.
