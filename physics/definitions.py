### A simple file to store constants and definitions

import numpy as np
import scipy
import matplotlib.pyplot as plt

"""Basic constants"""

h_planck = 4.135667696e-15 # eV s
cs = 29979245800 # cm /s 
rydberg_energy = 13.6056923 # eV
lyAlpha_energy = 10.198810993784 # eV
kB = 8.617343e-5 # eV / K
G = 6.674e-8 # dyne / cm^2 / g^2
pi = 3.14159265
e = 2.718281828

"""Define masses"""

mp = 1.67262171e-24 # grams
mp_eV = 938.272029e+6 # eV/c^2

me = 9.109383632e-28 # grams
me_eV = 510998.9461 # eV/c^2

"""Conversions"""

pc_to_cm = 3.0856776e+18 # cm
msol_to_grams = 1.9891e+33 # grams
eV_to_ergs = 1.60217653e-12 # ergs
GeV_To_InvSec = 1.52e+24 # sec
yrs_to_sec = 3.154e+7 

"""Cosmology""" 

# In acccordance to Planck 2018
h=0.674
H0=67.4/(pc_to_cm*10)
Omega_c = 0.12 / h**2
Omega_b = 0.0224 / h**2
Omega_m = Omega_c + Omega_b
Omega_lambda = 1 - Omega_m
Yp = 0.245

# Critical density
RhoCrit = 3*H0**2/(8*pi*G) * cs**2/eV_to_ergs # eV cm^-3

# Define Critical overdensity
Delta_c = 200

# Define Hydrogen number density
nH0 = RhoCrit*Omega_b*(1 - Yp)/mp_eV

# Define Helium number density
nHe0 = RhoCrit*Omega_b*(Yp)/mp_eV

# Define total number of atoms density
nA0 = nH0 + nHe0

## Define Hubble function

def Hubble(z):
    return H0*(Omega_lambda + Omega_m*(1+z)**3)**(1/2)

# Define IGM, CMB, and virial temperatures

def T_IGM(z):
    return 1/40 * (1+z)**2

def T_CMB(z):
    return 2.725*(1+z)

"""Define halo properties"""

def Tvir(z, mhalo):

    fac = (
            0.75*1.98e+4 * (1.22/0.6) *
            (Omega_m / Omega_t_m(z) * Delta_c / (18*pi**2) )**(1/3) *
            (mhalo/ 10**8)**(2/3) *
            (1+z)/10 * h**(2/3)
          )
    return fac

def Rvir(z, mhalo): 
    val = 323 * (mhalo/(1e+6))**(1/3) * ((1+z)/10)**(-1)
    return val

"""Define cosmological density parameters at time t"""

def Omega_t_m(z):
    return Omega_m*(1+z)**3 / (Hubble(z)/H0)**2

def Omega_t_lambda(z):
    return Omega_lambda / (Hubble(z)/H0)**2


"""Define growth factor"""

def growth_fac(z):
    
    # Define the constant denominator
    const = Omega_m*( Omega_m**(4/7) - Omega_lambda + (1 + Omega_m/2)*(1 + Omega_lambda/70) )**(-1)

    func = Omega_t_m(z)*( Omega_t_m(z)**(4/7) - Omega_t_lambda(z) + (1 + Omega_t_m(z)/2)*(1 + Omega_t_lambda(z)/70) )**(-1)

    return (func/const/(1+z))


"""
Collisional rate coefficients. units cm^3 / s
"""

# case B Hydrogen Recombination
def case_B(T):
    return 2.54e-13*(T/10**4)**(-0.8163)
    
# Electron attachment to H. Taken from Hirata 2006. Valid for T<=10^(4) K
def C_Hminus(T):
    return 3e-16*(T/300)**0.95 * e**(-T/9320)

# H2 formation via H minus
def C_H2(T):
    return 1.5e-9 * (T/300)**(-0.1)

"""
Define Lyman line energies
"""
def lyman_np_level(n):
    return rydberg_energy*(1 - 1/n**2)

"""
Define cross section
"""
def sigma_Hm(en):
    en = np.asanyarray(en)

    result = np.zeros_like(en, dtype=float)

    val = 7.928e+5 * h_planck**1.5
    
    mask = en>=0.755
    result[mask] = val*(en[mask] - 0.755)**1.5 / en[mask]**3
    
    return result

""" Define cooling rates """

def lambda_H2(xH2, nH, T): 
    """
    return H2 cooling rate: ergs cm^(-3) s^(-1)
    """
    return 3.41e-35 * T**(3.4) * xH2 * nH**2

## Define cooling rate according to Galli and Palla 1998
"""
This H2 cooling definition works in both low and LTE regimes. 
See examples.ipynb for a worked example and derivation. 
"""

def lambda_H2_cool (xH2, nH, T): 

    """
    Returns H2 cooling rate in units: ergs cm^-3 s^-1. 

    The function, however, is constructed first to calculate the cooling rate
    per H2 molecule. 

    First we define low density and LTE cooling rates and then bridge them together. 
    """

    # Define low density H cooling rate 
    def lambda_low_n(n,T):
        """ 
        Returns cooling rate per H2 molecule in low density regime. Depends on H number density.
        units: ergs s^(-1)
        """
        
        x=np.log10(T)
        val = 10**( - 103.0 + 97.59*x - 48.05*(x)**2 + 10.80 * (x)**3 - 0.9032*(x)**4 )
        return val * n 

    # Define LTE cooling rate
    def lambda_LTE(T): 

        """ 
        At LTE values, the cooling rate is independent of H number density 
        units: ergs s^(-1)
        """
        x = T/1e+3 
        # Rotational cooling term 
        lambda_rot_1 = ( (9.5e-22 * x**(3.76) ) / (1 + 0.12 * x**(2.1) ) *
                       np.exp(- (0.13/x)**3 ) 
                       )
        lambda_rot_2 = 3e-24*np.exp(-0.51 / x)

        # Vibrational Cooling Term 
        lambda_vib_1 = 6.7e-19*np.exp( - 5.86/x)
        lambda_vib_2 = 1.6e-18*np.exp( - 11.7/x)
    
        return (lambda_rot_1 + lambda_rot_2 + lambda_vib_1 + lambda_vib_2)

    # Bridge the two cooling rates together in ergs s^(-1)  
    lambda_H2_total = lambda_low_n (nH,T)*lambda_LTE(T)/(lambda_low_n (nH,T) + lambda_LTE(T))

    # Return total H2 cooling rate multiplied by H2 number density
    return lambda_H2_total * nH * xH2


# Define H cooling rate according to Cen 1992 

def lambda_H_cool(xe, nH, T): 
    val = 7.5e-19*(1 + (T/1e+5)**(1/2))**(-1) * np.exp(-118348/T)
    return val * xe * nH**2
