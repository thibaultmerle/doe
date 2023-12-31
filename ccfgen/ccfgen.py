#!/usr/bin/python3
'''
Simulation of Cross-Correlation Functions with an arbitrary number of components.
Each component is simulated by a Gaussian.
'''

import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt 

spectral_resolution = 60000 # Resolving power of the spectrograph
ncomp = 3                   # Number of radial velocity (RV) components in the CCF
vrads = [-40, -10, 80]      # Position of the RV components in km/s
ampls = [1.0, 0.7, 0.6]     # Normalized amplitudes of the RV components (the highest component being at 1.0) 
vsigs = [10, 10, 15]        # Intrinsic broadening of the RV components 
vmin, vmax = -400, 400	    # Velocity range in km/s
noise_offset = 0.0          
noise_amplitude = 0.1

ofn1 = f'ccf_{spectral_resolution}_{ncomp}_components.dat'
ofn2 = f'ccf_{spectral_resolution}_{ncomp}_components.npy'

velocity_elt = 299792.458/spectral_resolution
dv = round(velocity_elt/4, 1)

nrvs = round((vmax - vmin) /dv)
rvs = np.linspace(vmin, vmax, nrvs)

ccf = np.zeros(nrvs)
for i in range(ncomp):
    ccf +=  ampls[i] * np.sqrt(2*np.pi)*vsigs[i] * st.norm.pdf(rvs, vrads[i], vsigs[i])

noise = noise_offset + noise_amplitude*np.random.randn(nrvs)
noisy_ccf = ccf + noise

broad_kernel = st.norm.pdf(rvs, 0, velocity_elt)
convolved_noisy_ccf = np.convolve(noisy_ccf, broad_kernel, mode='same')
convolved_noisy_ccf /= max(convolved_noisy_ccf)

print(f'The spectral resolution is: {spectral_resolution}')
print(f'The velocity element is:    {velocity_elt:8.3f} km/s')
print(f'The velocity step is:       {dv:8.3f} km/s')

print(f'Velocity range: [{vmin}, {vmax}] km/s')
print(f'Number of velocity points: {nrvs}')

print(f'Output file1:  {ofn1}')
print(f'Output file2:  {ofn2}')
print(f'Control plot: {ofn1.replace(".dat", ".png")}')

data = np.array([rvs, convolved_noisy_ccf])
np.savetxt(ofn1, data.T, fmt='%12.4e')
np.save(ofn2, data)

plt.figure(1, figsize=(10,4), edgecolor='k', facecolor='w', tight_layout=True)
#plt.plot(rvs, noisy_ccf,'-+b')
plt.plot(rvs, ccf,'-k', label='Without instrumental broadening')
plt.plot(rvs, convolved_noisy_ccf,'-r', label='With instrumental broadening')
for i in range(ncomp):
	plt.axvline(x=vrads[i], ls='dotted', color='k', alpha=0.5, label=f'v{i+1} = {vrads[i]:.3f} km/s, $\\sigma$ = {vsigs[i]:.3f} km/s')
plt.xlabel('$v$ [km/s]')
plt.ylabel('CCF')
plt.legend(title=f'{ncomp} components')
plt.title(f'$R$ = {spectral_resolution}  Velocity element = {velocity_elt:.3f} km/s  Velocity sampling = {dv:.3f} km/s')
plt.xlim(vmin, vmax)
plt.savefig(ofn1.replace('.dat', '.png'))
plt.show()
