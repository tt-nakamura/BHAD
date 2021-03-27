from math import pi

G = 6.67384e-8 # gravitational const / cm^3/g/s^2
c = 2.99792458e10 # speed of light  / cm/s
hP = 6.6260688e-27 # Planck's const / erg s
hbar = hP/2/pi
kB = 1.3806488e-16 # Boltzmann const / erg/K
mH = 1.66053892e-24 # atomic mass unit / g
me = 9.1093898e-28 # electron mass / g
alphaEM = 1/137.0359895 # fine structure constant
sigmaT = 8*pi/3*(alphaEM*hbar/me/c)**2 # Thomson cross section / cm^2
a_rad =  pi**2/15*kB*(kB/hbar/c)**3 # radiation const / erg/cm^3/K^4
stefan = a_rad*c/4 # Stefan-Boltzmann const / erg/cm^2/s/K^4
eV = 1.60217733e-12 # 1 electron volt / erg
keV = 1e3*eV
Msun = 1.98892e33 # solar mass / g
LEsun = 4*pi*c*G*Msun*mH/sigmaT # Eddington luminosity / erg/s
MEsun = LEsun/c**2
pc = 3.08567758e18 # 1 parsec / cm
