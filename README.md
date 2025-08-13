# TAPPS
**T**hermophysical **A**b-initio **P**roperties from **P**hase **S**hifts

This repository contains a set of Python3 programs and corresponding data files enabling the calculation of thermophysical properties of noble gases from phase shift and bound state data.

A paper detailing the calculations is being finalized.

This dataset (v1.0.0) is now part of the [NIST Public Data Repository](https://data.nist.gov/od/id/mds2-3802).

## Main files

The main files are:

* ``thermophysicalPairProperties.py`` to compute the second virial coefficient, its first two temperature derivatives and the second acoustic virial coefficient for the pure and cross phases.
* ``Ne_avg_virials.py`` to compute the second virial coefficient, its first two temperature derivatives and the second acoustic virial coefficient for the normal mixture of neon isotopes.
* ``data/`` directory with compressed JSON files containing phase shifts and bound states computed for
  - <sup>3</sup>He, <sup>4</sup>He, <sup>3</sup>He–<sup>4</sup>He using the pair potential from [Czachorowski et al.](https://doi.org/10.1103/PhysRevA.102.042810)
  - <sup>20</sup>Ne, <sup>21</sup>Ne, <sup>22</sup>Ne and their cross values, using the pair potential from [Hellmann et al.](https://doi.org/10.1063/5.0047999)
  - <sup>40</sup>Ar using the pair potential from [Lang et al.](https://doi.org/10.1103/PhysRevA.109.052803)

## Dependencies

This software requires [NumPy](https://www.numpy.org) and [SciPy](https://www.scipy.org).   
Please see the [NumPy installation page](https://numpy.org/install/) or the [SciPy installation
page](https://scipy.org/install/) for details.

TAPPS has been developed using ``numpy 1.26`` and ``scipy 1.16``.

## Documentation

The ``README.TXT`` file contains detailed instructions on how to use the Python programs and the
classes contained therein. 

## Examples

* ``python3 thermophysicalPairProperties.py``: prints the the second virial coefficient, its first
  two temperature derivatives (multiplied by $T$ and $T^2$, respectively), and the second acoustic
  virial coefficient for <sup>4</sup>He at $T=273.16$ K. The output should be:
  ```
  B_He4(        273.16 )= 11.928088 cm3/mol
  TdBdT_He4(    273.16 )= -1.077099 cm3/mol
  T2d2BdT2_He4( 273.16 )= -0.748414 cm3/mol
  beta_a_He4(   273.16 )= 22.220466 cm3/mol
  ```
* ``python3 thermophysicalPairProperties.py -T 300``: prints the the second virial coefficient, its first
  two temperature derivatives (multiplied by $T$ and $T^2$, respectively), and the second acoustic
  virial coefficient for <sup>4</sup>He at $T=300$ K 
* ``python3 thermophysicalPairProperties.py -d data/He3_phase_shift_data.json.bz2``: prints the the
  second virial coefficient, its first two temperature derivatives (multiplied by $T$ and $T^2$,
  respectively), and the second acoustic virial coefficient for <sup>3</sup>He at $T=273.16$ K
* ``python3 thermophysicalPairProperties.py -d data/He3_phase_shift_data.json.bz2 -T 300``: prints the the
  second virial coefficient, its first two temperature derivatives (multiplied by $T$ and $T^2$,
  respectively), and the second acoustic virial coefficient for <sup>3</sup>He at $T=300$ K

More examples, including on how to compute the virial coefficients of the average neon mixture
and/or how to use TAPPS in your Python programs, are available in the ``README.TXT`` file.
