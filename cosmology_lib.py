import numpy as np
import parameter as par

#global Omega_M, Omega_L, Hubble

#Omega_M=0.3
#Omega_L=0.7
#Hubble=70.0

Omega_M0=par.Omega_M
#0.286000
#Omega_M0=0.25
Omega_L0=par.Omega_L
#Omega_L0=0.75
Hubble_h=par.Hubble/100.0
#Sigma_8=0.82
Sigma_8=0.8
Gamma_cosmo=0.21
Grav = 4.2994e-9 # [Mpc km^2/Msun/s^2]
H0=100.0

def z_hubbletime(t1):
    
    t= t1 * Hubble_h / 10.
    H0=1.02278  # H0=100km/s/Mpc
    a = (Omega_M0 / Omega_L0)**(1.0/3)*(np.sinh(3./2*np.sqrt(Omega_L0)*H0*t))**(2./3)
    z=1./a-1.
    
    return z

def func(z):
    Omega = Omega_L0 + Omega_M0
    e = np.sqrt(Omega_L0 + (1.0 - Omega) * (1+z)**2 + Omega_M0*(1+z)**3)
    
    f = 1./(1.+z)/e
    return f
    
def hubbletime_z(z):
    # approximate hubble time for a given redshift of a given cosmology
    H0=1.0 * 3.156/3.0857/10.
    
    if Omega_M0 <= 1.0:
        p1 = (1.-Omega_M0)/Omega_M0
        htime=2./(3.*np.sqrt(1.-Omega_M0)) * np.log(np.sqrt(p1) / (1.+z)**1.5 + np.sqrt(1.+p1/(1.+z)**3))
    else:
        htime = 2./3./(1.+z)/np.sqrt(1.+z)
    
    htime=htime/Hubble_h/H0
    return htime
    
def xH(z):
#;calculate hubble constant ( in physical units: km/s/Mpc ) at redshift z.
    z1 = 1.0 + z
    fac = Omega_L0 + (1.0 - Omega_L0 - Omega_M0)*z1*z1 + Omega_M0*z1*z1*z1
    xH = np.sqrt(fac);
    return xH

def Omega_m(z):
#// calculates the density parameter omega_m at redshift z. 
    omega = Omega_M0 * (1.0+z)*(1.0+z)*(1.0+z) / xH(z)**2
    return omega

def rho_crit(z):
#// critical density of the universe
    h = xH(z) * H0
    rho_c = 3*h*h / (8*np.pi*Grav) #   //in unit of Msun h^-1/(Mpc h^-1)^3
    return rho_c


def delta_vir(z):
#// calculate the virial density in terms of critical density of the Universe.
#// We use the fitting formulae of Bryan & Norman, 1998, ApJ, 495, 80.
#// These are accurate to better than 1% for 0.1<Omega_M0<1.0
    dc = 200.
    o_m = Omega_m(z)
    x= o_m - 1.
    if (Omega_M0 > 1. and Omega_L0 == 0.0):
        dc=18.0*np.pi*np.pi + 60.0*x - 32.0*x*x
    if (np.abs(1.0-Omega_M0-Omega_L0) <  1.0E-4):
        dc=18.0*np.pi*np.pi + 82.0*x - 39.0*x*x
    return dc


