import numpy as np
from scipy.special import expn
from constants import c,hP,kB

kappa_es = 0.4 # opacity by electron scattering / cm^2/g
kappa_0 = 6.45e22 # krammers opacity at (rho,T)=(1,1) in cgs

def kappa_ff(rho,T):
    """ krammers law for opacity by free-free absorption
    rho = density / g/cm^3
    T = temperature / K
    return opacity / cm^2/g
    """
    return kappa_0*rho/T**3.5

nu_0 = 2.4e21*15 # electron-ion coupling (Coulomb log = 15)

def nu_E(rho,T):
    """ collision rate of ions with electrons 
    rho = density / g/cm^3
    T = electron temperature / K
    return collision rate / s^-1
    """
    return nu_0*rho/T**1.5

def BlackBody(nu,T):
    """ black-body spectrum
    nu = frequency / Hz
    T = temperature / K
    return d(intensity)/dln(nu) / erg/s/cm^2/steradian
    """
    return 2*hP*nu**4/c**2/np.expm1(hP*nu/kB/T)

def InvCompton(nu,T,Es):
    """ spectrum of inverse-compton scattered photons
    assuming compton y parameter is unity
    nu = frequency / Hz
    T = temperature of electrons / K
    Es = energy of photons before scattering / erg
    return: normalized spectrum
      dI/dln(nu) /[\int dI/d(nu) d(nu)] if nu >= Es/hP
      0 if nu < Es/hP
    reference: S.L.Shapiro, A.P.Lightman and D.M.Eardley
    The Astrophysical Journal 204 (1976) 187, eq.(17)
    """
    x = hP*nu/kB/T
    y = Es/kB/T
    I = (1 + x*(1 + x/2*(1 + x/3*(1 + x/4))))*np.exp(-x)
    I /= (50 + y*(26 + y*(7 + y)))/24*np.exp(-y) + expn(1,y)
    return np.where(x<y, 0, I)
