#!/usr/bin/env python3

#
# thermophysicalPairProperties.py
#
# A program to calculate thermophysical properties from phase shifts
# (c) Giovanni Garberoglio, 2025
# garberoglio@ectstar.eu
#

#
# If useful, please cite
# https://doi.org/DOI/OF/OUR/PAPER
#

import numpy as np
import scipy as sp
import json
import bz2
import os
import sys

kB   = sp.constants.k
NA   = sp.constants.N_A
hbar = sp.constants.hbar
amu  = sp.constants.physical_constants['atomic mass constant'][0]

def load_json(filename):
    # Check if the file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File '{filename}' does not exist.")

    # Try to load the file content as JSON
    try:
        if filename.endswith('.bz2'):
            # Open with bz2 if compressed
            with bz2.open(filename, 'rt') as file:
                data = json.load(file)
        else:
            # Otherwise, open as a regular JSON file
            with open(filename, 'r') as file:
                data = json.load(file)
    except json.JSONDecodeError as e:
        raise ValueError(f"File '{filename}' is not a valid JSON file.") from e

    return data


class ThermophysicalPairProperties():

    def __init__(self, filename, verbose=False):
        data_dict = load_json(filename)
        self.data_dict = data_dict
        self.verbose = verbose
        
        mu = data_dict['reduced mass']['value']
        self.mu = mu
        self.xi = hbar*hbar/(2.0*mu*amu*1e-20*kB)

        # min and max temperatures for which virial coefficients are reliable
        self.Tmin = data_dict['Tmin']
        self.Tmax = data_dict['Tmax']        
        
        if data_dict['identical particles'] == True:
            I = data_dict['nuclear spin']
            self.I = I
        else:
            self.I = -1
            
        bound_data = np.array(data_dict['bound states']['data'])
        if len(bound_data) > 0:
            self.L_bound = bound_data[:,0]
            self.E_bound = bound_data[:,1]
            if self.I >= 0:                
                self.f_bound = np.power(-1,self.L_bound+2*I)/(2*I+1)
            else:
                self.f_bound= np.zeros_like(self.L_bound)                
        else:
            self.L_bound = []
            self.E_bound = []
            self.f_bound = 0
            
        self.Eval = np.array(data_dict['energies']['data'])
        self.kval = np.sqrt(self.Eval/self.xi)
        self.L = np.array(data_dict['angular momenta']['data']) 
        
        if self.I >= 0:            
            f = np.power(-1,self.L+2*I)/(2*I+1)
            self.f = np.power(-1,2*I)/(2*I+1)
        else:
            f = np.zeros_like(self.L)
            self.f  = 0.0
        self.d_L = (1+f)    
        
        self.phase_shifts = np.array(data_dict['phase shifts']['data'])

        self.integration='Gaussian'
        # gk_n is the number of points for fixed Gaussian integration
        # 3000 seems already to be enough for virial and transport coefficients, so here we try to
        # make it future proof (just in case)
        self.gk_n = 5000
        # following are parameters for adaptive integration. However, it does not perform as well as the
        # Gaussian integration above, especially for transport coefficients, at least when using
        # scipy 1.15.2
        # hence we default to Gaussian integration
        self.epsrel = 1e-7
        self.limit  = 512
        
    def spline_integrate(self, x, y):
        # \int y(x) dx using cubic-spline representation of y
        ys = sp.interpolate.CubicSpline(x,y)
        q = sp.integrate.quad(lambda t: ys(t), x[0],x[-1],
                                    epsrel=self.epsrel, limit=self.limit)
        return(q[0])

    def spline_integrate_gk(self, x, y):
        ys = sp.interpolate.CubicSpline(x,y)
        q = sp.integrate.fixed_quad(lambda t: ys(t), x[0],x[-1],
                                    n=self.gk_n)
        return(q[0])
        
    def In(self, n, T):
        # returns \int S(k) k^{n+1} \exp(-xi k^2/T) dk
        # which is NOT the I_n of the paper
        # I_n(paper) = -N Lambda^5/(2 pi)^2 * In
        
        d = (2*self.L+1) * self.d_L
        d = d.reshape(-1,1)
        
        kbf = np.power(self.kval,n+1) * np.exp(-self.Eval/T)
        In_integrand = np.sum(d * self.phase_shifts * kbf,axis=0)

        # We know that the integrand is 0 and K=0
        K = np.insert(np.copy(self.kval), 0, 0.0)
        In_integrand = np.insert(In_integrand, 0, 0.0)

        if self.verbose == True:
            print('K',K[0:5])
            print('In',In_integrand[0:5])
        
        match self.integration:
            case 'Gaussian':
                In = self.spline_integrate_gk(K, In_integrand)
            case _:
                In = self.spline_integrate(K, In_integrand)
                
        return(In)

    def check_temperature(self,T):
        if T < self.Tmin or T > self.Tmax:
            species_name = self.data_dict['metadata']['species']
            errormsg = "Error: Can compute " + species_name + " virial for "
            errormsg += str(self.Tmin)+" <= T <= "+str(self.Tmax)+"\n"
            errormsg += str(T)+" is out of bounds"
                    
            raise ValueError(errormsg)

    def second_virial(self, T):
        """
        Computes the second virial coefficients from phase shifts.
        Input: T in K
        Output: B(T) in cm3/mol
        """
        self.check_temperature(T)
        A = NA * 1e-24 # factor for converting to cm3/mol
        
        # de Broglie thermal wavelength
        Lambda = np.sqrt(self.xi * 4.0 * np.pi/T)
        Lambda3 = np.power(Lambda, 3.0)

        # thermal contribution
        I0 = self.In(0, T)        
        B_th  = I0
        B_th *= -A*Lambda3*(Lambda/(2*np.pi))**2
        
        # ideal gas contribution
        B_ideal = -A*self.f*Lambda3/16.0
    
        # bound states contribution
        if len(self.E_bound) > 0.0:
            L = self.L_bound
            E = self.E_bound
            d = (2*L+1)*(1+self.f_bound)
            B_bound = np.sum(d * np.expm1(-E/T))
            B_bound *= -A*0.5*Lambda3
        else:
            B_bound = 0.0

        if self.verbose == True:
            print('B_th(%g) = %g' % (T,B_th))
            print('B_ideal(%g) = %g' % (T,B_ideal))
            print('B_bound(%g) = %g' % (T,B_bound))

        return(B_th + B_ideal + B_bound)

    def B(self,T):
        return(self.second_virial(T))

    def TdBdT(self, T):
        """
        Computes T times the temperature derivative of the second virial coefficient from phase shifts.
        Input: T in K
        Output: T dB(T)/dT in cm3/mol
        """
        self.check_temperature(T)        
        A = NA * 1e-24 # factor for converting to cm3/mol
        
        # de Broglie thermal wavelength
        Lambda = np.sqrt(self.xi * 4.0 * np.pi/T)
        Lambda3 = np.power(Lambda, 3.0)

        # thermal contribution
        I0 = self.In(0, T)    
        I2 = self.In(2, T)
        TdBdT_th = -2.5*I0 + I2*Lambda*Lambda/(4*np.pi)
        TdBdT_th *= -A*Lambda3*(Lambda/(2*np.pi))**2
        
        # ideal gas contribution
        B_ideal = -A*self.f*Lambda3/16.0
        TdBdT_ideal = -3*B_ideal/2
        
        # bound-state contribution
        if len(self.E_bound) > 0.0:
            L = self.L_bound
            E = self.E_bound
            d = (2*L+1)*(1+self.f_bound)
            x = E/T
            eem1 = np.expm1(-x)
            ee = eem1 + 1
            TdBdT_bound = np.sum(d * (-1.5*eem1 + x*ee))
            TdBdT_bound *= -A*0.5*Lambda3
        else:
            TdBdT_bound = 0.0

        if self.verbose == True:
            print('TdBdT_th(%g) = %g' % (T,TdBdT_th))
            print('TdBdT_ideal(%g) = %g' % (T,TdBdT_ideal))
            print('TdBdT_bound(%g) = %g' % (T,TdBdT_bound))

        return(TdBdT_th + TdBdT_ideal + TdBdT_bound)

    def T2d2BdT2(self, T):
        """
        Computes T times the temperature derivative of the second virial coefficient from phase shifts.
        Input: T in K
        Output: T dB(T)/dT in cm3/mol
        """
        self.check_temperature(T)        
        A = NA * 1e-24 # factor for converting to cm3/mol
        
        # de Broglie thermal wavelength
        Lambda = np.sqrt(self.xi * 4.0 * np.pi/T)
        Lambda3 = np.power(Lambda, 3.0)

        # thermal contribution
        I0 = self.In(0, T)    
        I2 = self.In(2, T)
        I4 = self.In(4, T)
        T2d2BdT2_th  = 35*I0/4
        T2d2BdT2_th -= 7*I2*Lambda*Lambda/(4*np.pi)
        T2d2BdT2_th += I4*(Lambda*Lambda/(4*np.pi))**2
        T2d2BdT2_th *= -A*Lambda3*(Lambda/(2*np.pi))**2
        
        # ideal gas contribution
        B_ideal = -A*self.f*Lambda3/16.0
        T2d2BdT2_ideal = 15*B_ideal/4

        # bound-state contribution
        if len(self.E_bound) > 0.0:
            L = self.L_bound
            E = self.E_bound
            d = (2*L+1)*(1+self.f_bound)
            x = E/T
            eem1 = np.expm1(-x)
            ee = eem1 + 1
            T2d2BdT2_bound = np.sum(d * ( 15*eem1/4.0 + (-5*x+x*x)*ee ))
            T2d2BdT2_bound *= -A*0.5*Lambda3
        else:
            T2d2BdT2_bound = 0.0

        if self.verbose == True:
            print('T2d2BdT2_th(%g) = %g' % (T,T2d2BdT2_th))
            print('T2d2BdT2_ideal(%g) = %g' % (T,T2d2BdT2_ideal))
            print('T2d2BdT2_bound(%g) = %g' % (T,T2d2BdT2_bound))

        return(T2d2BdT2_th + T2d2BdT2_ideal + T2d2BdT2_bound)

    def second_acoustic_virial(self, T):
        """
        Computes the second acoustic virial coefficient
        """
        gamma0 = 5.0/3.0
        B = self.B(T)
        TdBdT = self.TdBdT(T)
        T2d2BdT2 = self.T2d2BdT2(T)

        beta_a  = 2.0*B
        beta_a += 2.0*(gamma0-1.0)*TdBdT;
        beta_a += (gamma0-1.0)**2 * T2d2BdT2/gamma0

        return(beta_a)

    def beta_a(self, T):
        return(self.second_acoustic_virial(T))

    
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Compute thermophysical pair properties from phase shifts')
    parser.add_argument('-d', default='data/He4_phase_shift_data.json.bz2', help='phase shift data file')
    parser.add_argument('-T', default='273.16', help='temperature[s] (number | python list | numpy expression)')
    parser.add_argument('-v', type=bool, default=False, help='be verbose')
    
    args=parser.parse_args()
    datafile=args.d
    verbose=args.v

    T=eval(args.T)
    if type(T)==int or type(T)==float:
        T = np.array([T])
        

    if len(datafile) == 0:
        print("Need a data file (-d option)")
        exit(3)
        
    pp = ThermophysicalPairProperties(datafile, verbose=verbose)
    species = pp.data_dict['metadata']['species']

    for t in T:
        print("B_%s(        %g )= %.6f cm3/mol" % (species, t, pp.B(t)))
        print("TdBdT_%s(    %g )= %.6f cm3/mol" % (species, t, pp.TdBdT(t)))
        print("T2d2BdT2_%s( %g )= %.6f cm3/mol" % (species, t, pp.T2d2BdT2(t)))
        print("beta_a_%s(   %g )= %.6f cm3/mol" % (species, t, pp.beta_a(t)))        
        print("")
