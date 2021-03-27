import numpy as np
import matplotlib.pyplot as plt
from BHAD import SLEModel

M = 10 # BH mass / solar mass
M_dot = 1 # accretion rate / critical rate
alpha = 1
r_max = 3e3 # radius of disk's outer edge / 2GM/c^2

plt.figure(figsize=(6,5))
plt.subplots_adjust(top=0.98, bottom=0.1,
                    right=0.98, left=0.12, hspace=0)

s = SLEModel(M, M_dot, alpha, r_max)
r = s.radius()
Ti = s.ion_temperature()
Te = s.electron_temperature()
Teff = s.effective_temperature()
plt.subplot(2,1,1)
plt.loglog(r, Ti, 'b', label='ion')
plt.loglog(r, Te, 'r--', label='electron')
plt.loglog(r, Teff, 'g', label='effective')
plt.ylabel('temperature / K', fontsize=14)
plt.tick_params('x', which='both', direction='in')
plt.xlim([s.r_in, s.r_out])
plt.legend()

p = s.pressure()
p_gas = s.gas_pressure()
plt.subplot(2,1,2)
plt.loglog(r, p, 'b')
plt.loglog(r, p_gas, 'r--', label='gas pressure')
plt.ylabel('pressure / dyne/cm$^2$', fontsize=14)
plt.xlim([s.r_in, s.r_out])
plt.legend()

plt.xlabel('$r$ = distance from black hole / cm', fontsize=14)
plt.savefig('fig4.eps')
plt.show()
