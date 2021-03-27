import numpy as np
import matplotlib.pyplot as plt
from BHAD import StandardModel
from radiation import kappa_es

M = [10]*4 # BH mass / solar mass
M_dot = [1, 0.1, 1, 0.1] # accretion rate / critical rate
alpha = [1, 1, 0.1, 0.1]
r_max = 3e3 # radius of disk's outer edge / 2GM/c^2

color = ['g', 'orange', 'b', 'r']

plt.figure(figsize=(6,8))
plt.subplots_adjust(top=0.98, bottom=0.1,
                    right=0.98, left=0.12, hspace=0)

for M,M_dot,alpha,c in zip(M,M_dot,alpha,color):
    s = StandardModel(M, M_dot, alpha, r_max)
    r = s.radius()
    T_eff = s.effective_temperature()
    p = s.pressure()
    p_gas = s.gas_pressure()
    kappa = s.opacity()
    label = r'$\dot m=%s, \alpha=%s$'%(str(M_dot),str(alpha))
    plt.subplot(3,1,1);
    plt.loglog(r, T_eff, c=c, label=label)
    plt.subplot(3,1,2);
    plt.loglog(r, p, c=c)
    l = plt.loglog(r, p_gas, c=c, ls='--')
    if c=='b': l[0].set_label('gas pressure')
    plt.subplot(3,1,3);
    plt.loglog(r, kappa, c=c)
    l = plt.loglog(r, [kappa_es]*len(r), c=c, ls='--')
    if c=='r': l[0].set_label('electron scattering')

y_label = ['$T_{eff}$ / K',
           'pressure / dyne/cm$^2$', 'opacity / cm$^2$/g']
y_ticks = [[1e5,1e6], [1e9,1e12,1e15], [1e-1,1e0]]

for i in range(3):
    plt.subplot(3,1,i+1)
    plt.ylabel(y_label[i], fontsize=14)
    plt.yticks(y_ticks[i])
    if i<2: plt.tick_params('x', which='both', direction='in')
    plt.xlim([s.r_in, s.r_out])
    plt.legend()

plt.xlabel('$r$ = distance from black hole / cm', fontsize=14)
plt.savefig('fig2.eps')
plt.show()
