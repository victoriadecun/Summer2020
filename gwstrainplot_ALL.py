import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from collections import OrderedDict
#%matplotlib inline                                                                                   
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)

#local imports                                                                                        
from detection_rate_main import detection_rate_main, z_at
from utils.mergerrate import MergerRate
from utils.mbhbinaries import  mass_ratio_func
from utils.evolveFDFA import EvolveFDFA
from utils.basicmergers import NoDelay, ConstantTime, NoDelayOrigExtract

import sys
sys.path.append('utils/evolve_lzk/zcode/')
#from utils.evolve_lzk.evolveLZK import EvolveLZK                                                     

from gwsnrcalc.utils.waveforms import PhenomDWaveforms
from gwsnrcalc.gw_snr_calculator import snr
from gwsnrcalc.utils.readnoisecurves import read_noise_curve
import scipy.constants as ct
from scipy import stats
Msun=1.989e30
import pdb
import os

#GW strain plot                                                                                       

fig, ax = plt.subplots(1,1)
fig.set_size_inches(10,8)
phenomdwave = PhenomDWaveforms(num_points=4096)

q=0.2
z=1.0
s=0.8
start_time = 100.0
end_time = 0.0
for M, z, q, letter in [[5.64555e6, 1.9163191, 0.0089485200, 'A'], [7.85521e6, 0.39280810, 0.25782200 'B'], [1.48301e6, 10.762340, 0.19103100, 'C'], [4.66498e6, 9.4247567, 0.29180500, 'D'], [82245.4, 3.7944614, 0.47851300, 'E'], [8.29909e6, 1.7536451, 0.67958900, 'F'], [676377, 11.378336, 0.86115400, 'G'], [2.58093e6, 9.6918367, 0.44843200, 'H'], [2.24031e6, 4.6032768, 0.30951800, 'I'], [1.70916e6, 8.9868174, 0.34992000, 'J'], [336181.25, 5.2838543, 0.025792000, 'K']]:
    ax.set_ylim(6e-23, 5e-14)
    ax.set_xlim(1e-6, 1e-3)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(18)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(18)

    m1 = M/(1+q)
    m2 = M*q/(1+q)
    wave = phenomdwave(m1, m2, s, s, z, start_time, end_time)

    f = wave.freqs
    hc = wave.hc*np.sqrt(96/200) #averaging factor see Cornish and Robson 2018                         

    ins = np.where(f<wave.fmrg)[0]
    mrg = np.where((f>=wave.fmrg)&(f<=wave.fpeak))[0]
    rd = np.where(f>wave.fpeak)[0]

    ax.loglog(f[ins], hc[ins], color='blue', lw = 2)
    ax.loglog(f[mrg], hc[mrg], color='green', lw = 2)
    ax.loglog(f[rd], hc[rd], color='brown', lw = 2)

    ax.text(f[rd[0]]*1.1, hc[rd[0]]*1.1, letter, fontsize=18)

    N = m1*m2/(m1+m2)**2.
    start_times = np.array([100.0, 10.0,1.0])*ct.Julian_year/(1+z)

    tau = N*(start_times*ct.c)/(5.*(m1+m2)*Msun*ct.G/(ct.c**2.))
    flow = 1./(8.*ct.pi*(m1+m2)*Msun*ct.G/(ct.c**2.)*tau**(3./8.))*(1.+((11./32)*N+743./2688.)*tau**(-
1./4.))*ct.c/(1+z)

    vert = [45, 30, 15]
    for i,(flow_i, st) in enumerate(zip(flow, start_times)):
        f_ind = np.where(f >= flow_i)[0][0]
        val = st*(1+z)/ct.Julian_year
        string = "%i yrs"%int(val) if val != 1.0 else "%i yr"%int(val)
        ax.annotate(string,xy=(f[f_ind], hc[f_ind]), xycoords='data',xytext=(0,vert[i]), ha='center', 
textcoords='offset points',arrowprops=dict(arrowstyle="-",linewidth = 1.75), fontsize=14)

ax.plot(1,1, color='blue', ls='solid', label='Inspiral')
ax.plot(1,1, color='green', ls='solid', label='Merger')
ax.plot(1,1, color='brown', ls='solid', label='Ringdown')
for nc in ['PL', 'HB_wd_noise']:
    fn, hn = read_noise_curve(nc, noise_type_out='char_strain')
    if nc == 'HB_wd_noise':
        nc = 'Galactic Background'
    ax.loglog(fn, hn, label=nc, ls='dashed')

ax.tick_params(labelsize=16)
ax.set_xlabel('Frequency (Hz)', fontsize=20)
ax.set_ylabel(r'Characteristic Strain', fontsize=20)
ax.set_xlim(2e-7, 1e0)
ax.legend(prop={'size':14})
fig.savefig('gwstrainplot_ALL.pdf', dpi=200)
