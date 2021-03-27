import numpy as np
import matplotlib.pyplot as plt
from BHAD import StandardModel

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
    H = s.height()
    u = s.infall_velocity()
    v = s.rotation_velocity()
    tau = s.optical_depth()
    label = r'$\dot m=%s, \alpha=%s$'%(str(M_dot),str(alpha))
    plt.subplot(3,1,1); plt.loglog(r, u/v, c=c)
    plt.subplot(3,1,2); plt.loglog(r, tau, c=c)
    plt.subplot(3,1,3); plt.loglog(r, H/r, c=c, label=label)

y_label = ['$u/v$', 'optical depth', '$H/r$']
y_ticks = [[1e-5,1e-4,1e-3], [1e2,1e3], [1e-3,1e-2]]

for i in range(3):
    plt.subplot(3,1,i+1)
    plt.ylabel(y_label[i], fontsize=14)
    plt.yticks(y_ticks[i])
    if i<2: plt.tick_params('x', which='both', direction='in')
    plt.xlim([s.r_in, s.r_out])

plt.xlabel('$r$ = distance from black hole / cm', fontsize=14)
plt.legend()
plt.savefig('fig1.eps')
plt.show()
