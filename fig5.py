import numpy as np
import matplotlib.pyplot as plt
from BHAD import StandardModel,SLEModel
from constants import keV

M = 10 # BH mass / solar mass
M_dot = 1 # accretion rate / critical rate
alpha = 1
r_out = 3e3 # radius of disk's outer edge / 2GM/c^2
nu = np.geomspace(1e13, 1e21, 128) # frequency / Hz

s = StandardModel(M, M_dot, alpha, r_out)
t = SLEModel(M, M_dot, alpha, r_out)
plt.loglog(nu, s.spectrum(nu), 'g', label='standard model')
plt.loglog(nu, t.spectrum(nu, 0.5*keV), 'r--', label='hot model')
plt.axis([nu[0], nu[-1], 2e-11, 2e-1])
plt.xlabel(r'$\nu$ = frequency / Hz', fontsize=14)
plt.ylabel(r'$\nu S_\nu$ = brightness / erg/s/cm$^2$', fontsize=14)
plt.legend()
plt.tight_layout()
plt.savefig('fig5.eps')
plt.show()
