#!/usr/bin/env python3

#
# Ne_avg_virials.py
#
# A program to calculate thermophysical properties of the average Neon mixture
# (c) Giovanni Garberoglio, 2024
# garberoglio@ectstar.eu
#
# If useful, please cite
# https://doi.org/DOI/OF/OUR/PAPER
#

import argparse
import itertools as it
import numpy as np
import thermophysicalPairProperties as tpp

Ne_mol_fractions = {'Ne20':0.9048 ,
                    'Ne21':0.0027,
                    'Ne22':0.0925}

class Ne_Average_Mixture():

    def __init__(self, verbose=False):
        # prepare the weights of the pure and cross virials
        self.yij = []
        identifiers = []

        for s1,s2 in it.combinations_with_replacement(Ne_mol_fractions.keys(),2):
            y1 = Ne_mol_fractions[s1]
            y2 = Ne_mol_fractions[s2]
            if s1 == s2:
                self.yij += [ y1*y2 ]
                identifiers += [ s1 ]
            else:
                self.yij += [ 2*y1*y2 ]
                identifiers += [ s1+"-"+s2 ]

        if verbose == True:
            for w,x in zip(yij, identifiers):
                print("weight of B_%9s = %g" % (x, w))
                print("")
    
        suffix = "_phase_shift_data.json.bz2"
        self.pp     = [tpp.ThermophysicalPairProperties(x+suffix) for x in identifiers]
        
    def B(self,T):
        B_data = [x.B(T) for x in self.pp]
        B      = np.dot(self.yij, B_data)
        return(B)
        
    def TdBdT(self, T):
        TdBdT_data    = [x.TdBdT(T) for x in self.pp]
        TdBdT    = np.dot(self.yij, TdBdT_data)
        return(TdBdT)

    def T2d2BdT2(self, T):
        T2d2BdT2_data = [x.T2d2BdT2(T) for x in self.pp]
        T2d2BdT2 = np.dot(self.yij, T2d2BdT2_data)    
        return(T2d2BdT2)

    def beta_a(self, T):
        beta_a_data   = [x.beta_a(T) for x in self.pp]
        beta_a   = np.dot(self.yij, beta_a_data)
        return(beta_a)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute thermophysical pair properties of the average Ne mixture from phase shifts')
    parser.add_argument('-T', default='273.16', help='temperature[s] (number | python list | numpy expression)')
    parser.add_argument('-v', type=bool, default=False, help='be verbose')
    
    args=parser.parse_args()
    verbose=args.v

    T=eval(args.T)
    if type(T)==int or type(T)==float:
        T = np.array([T])

    NeAvg = Ne_Average_Mixture(verbose)

    for t in T:
        B        = NeAvg.B(t)
        TdBdT    = NeAvg.TdBdT(t)
        T2d2BdT2 = NeAvg.T2d2BdT2(t)   
        beta_a   = NeAvg.beta_a(t)
    
        print("B_Neavg(        %g )= %.6f cm3/mol" % (t, B))
        print("TdBdT_Neavg(    %g )= %.6f cm3/mol" % (t, TdBdT))
        print("T2d2BdT2_Neavg( %g )= %.6f cm3/mol" % (t, T2d2BdT2))
        print("beta_a_Neavg(   %g )= %.6f cm3/mol" % (t, beta_a))        
        print("")
    
