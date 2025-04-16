thermophysicalPairProperties.py is a Python3 program and class to compute the second virial
coefficient B(T), its first and second temperature derivative, and the second acoustic virial
coefficient of atomic gases or mixtures from phase-shift data.

Description of the files:

*_phase_shift_data.json.bz2:
    phase shift and bound state data, in compressed JSON format for
    Pure gases: He3, He4, Ne20, Ne21, Ne22, Ar40
    Mixtures: He3-He4, Ne20-Ne21, Ne20-Ne22, Ne21-Ne22

	The phase shifts are enough to provide values of the virial coefficients within the published
	uncertainties up to a temperature Tmax = 1500 K and beyond a minimum temperature of
	 0.01 K for He3, He4, and their mixture
	 1 K    for Ne20, Ne21, Ne22, and their mixtures
	10 K    for Ar40
    A python ValueError will be raised if the requested temperature is outside these limits.    
    
thermophysicalPairProperties.py:
    main python3 program and class. Needs at least one of the data files above.
    See "USAGE AS CLASS (in your python programs)" below for instructions
    
Ne_avg_virials.py:
    main python3 program and class to compute the virial coefficients for the normal neon
    mixture. Needs all the Ne*bz2 data files and
    thermophysicalPairProperties.py
    See "USAGE AS CLASS | Neon mixture" below for instructions
    
--------------------------------------------------------------------------------

USAGE AS PROGRAM
----------------

python3 thermophysicalPairProperties.py [-T temperature -d datafile -h ]
* 'temperature' is the required temperature
(advanced usage: this can be any valid python iterable returning temperatures, e.g. -T "[1,10,100,1000]")
* 'datafile' is one of the data files above
* '-h' prints help on the available options

python3 Ne_avg_virials.py [-T temperature -h ]
* 'temperature' is the required temperature
(advanced usage: this can be any valid python iterable returning temperatures, e.g. -T "[1,10,100,1000]")
* '-h' prints help on the available options

EXAMPLES
--------

$ denotes the input prompt. The subsequent lines are the output.

$ python3 thermophysicalPairProperties.py
B_He4(        273.16 )= 11.928088 cm3/mol
TdBdT_He4(    273.16 )= -1.077099 cm3/mol
T2d2BdT2_He4( 273.16 )= -0.748414 cm3/mol
beta_a_He4(   273.16 )= 22.220466 cm3/mol

$ python3 thermophysicalPairProperties.py -d data/He3_phase_shift_data.json.bz2 
B_He3(        273.16 )= 12.028842 cm3/mol
TdBdT_He3(    273.16 )= -1.205720 cm3/mol
T2d2BdT2_He3( 273.16 )= -0.455542 cm3/mol
beta_a_He3(   273.16 )= 22.328580 cm3/mol

$ python3 thermophysicalPairProperties.py -T '[1,10,100,1000]'
B_He4(        1 )= -475.697195 cm3/mol
TdBdT_He4(    1 )= 670.284965 cm3/mol
T2d2BdT2_He4( 1 )= -1793.712314 cm3/mol
beta_a_He4(   1 )= -536.004386 cm3/mol

B_He4(        10 )= -23.123061 cm3/mol
TdBdT_He4(    10 )= 41.020085 cm3/mol
T2d2BdT2_He4( 10 )= -82.473578 cm3/mol
beta_a_He4(   10 )= -13.545630 cm3/mol

B_He4(        100 )= 11.675048 cm3/mol
TdBdT_He4(    100 )= 2.090648 cm3/mol
T2d2BdT2_He4( 100 )= -6.998566 cm3/mol
beta_a_He4(   100 )= 24.271343 cm3/mol

B_He4(        1000 )= 9.550472 cm3/mol
TdBdT_He4(    1000 )= -2.255590 cm3/mol
T2d2BdT2_He4( 1000 )= 1.992370 cm3/mol
beta_a_He4(   1000 )= 16.624790 cm3/mol

$ python3 Ne_avg_virials.py 
B_Neavg(        273.16 )= 11.018930 cm3/mol
TdBdT_Neavg(    273.16 )= 6.865094 cm3/mol
T2d2BdT2_Neavg( 273.16 )= -17.730464 cm3/mol
beta_a_Neavg(   273.16 )= 26.463194 cm3/mol

$ python3 Ne_avg_virials.py -T '[10,100,1000]'
B_Neavg(        10 )= -909.361067 cm3/mol
TdBdT_Neavg(    10 )= 2532.397550 cm3/mol
T2d2BdT2_Neavg( 10 )= -11578.183779 cm3/mol
beta_a_Neavg(   10 )= -1529.707743 cm3/mol

B_Neavg(        100 )= -4.288236 cm3/mol
TdBdT_Neavg(    100 )= 27.746620 cm3/mol
T2d2BdT2_Neavg( 100 )= -64.375995 cm3/mol
beta_a_Neavg(   100 )= 11.252090 cm3/mol

B_Neavg(        1000 )= 13.913426 cm3/mol
TdBdT_Neavg(    1000 )= -0.692441 cm3/mol
T2d2BdT2_Neavg( 1000 )= -1.959108 cm3/mol
beta_a_Neavg(   1000 )= 26.381170 cm3/mol

USAGE AS A CLASS (in your python programs)
------------------------------------------

import thermophysicalPairProperties as tpp

datafile='data/He4_phase_shift_data.json.bz2'
pp = tpp.ThermophysicalPairProperties(datafile)

# set the temperature
T = 300
# compute properties from the phase shifts
B      = pp.B(T)        # second virial coefficient in cm3/mol
TB     = pp.TdBdT(T)    # first derivative multiplied by T in cm3/mol
T2B    = pp.T2d2BdT2(T) # second derivative multiplied by T^2 in cm3/mol
beta_a = pp.beta_a(T)   # second acoustic virial coefficient in cm3/mol

print("B_He4(%g)      = %g cm3/mol" % (T, B))
print("beta_a_He4(%g) = %g cm3/mol" % (T, beta_a))

USAGE AS A CLASS | Neon mixture
-------------------------------

import Ne_avg_virials as na

# instantiate the class
NeAvg = na.Ne_Average_Mixture()

T=300
# compute the properties
Bmix = NeAvg.B(T)

print("Neon mixture B(%g) = %g cm3/mol" % (T, Bmix)) 
