import numpy as np
import matplotlib.pyplot as plt
from BHAD import StandardModel

M = [10,10,100,100] # BH mass / solar mass
M_dot = [1, 0.1, 1, 0.1] # accreation rate / M_dot_crit
alpha = 0.1
r_max = 3e3 # radius of disk's outer edge / 2GM/c^2
nu = np.geomspace(1e13, 1e18, 256)

color = ['g', 'orange', 'b', 'r']

for M,M_dot,c in zip(M,M_dot,color):
    s = StandardModel(M, M_dot, alpha, r_max)
    I = s.spectrum(nu)
    label = r'$M=%s, \dot m=%s$'%(str(M),str(M_dot))
    plt.loglog(nu, I, c=c, label=label)

plt.ylabel(r'$\nu S_\nu$ = brightness / erg/s/cm$^2$', fontsize=14)
plt.axis([nu[0], nu[-1], 2e-11, 2e-1])
plt.legend()

plt.xlabel(r'$\nu$ = frequcency / Hz', fontsize=14)
plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
