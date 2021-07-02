import numpy as np
import matplotlib.pyplot as plt


from evolveFDFA import EvolveFDFA
from gwsnrcalc.utils.waveforms import PhenomDWaveforms
from gwsnrcalc.gw_snr_calculator import snr
from gwsnrcalc.utils.readnoisecurves import read_noise_curve

#load data
#data = np.genfromtxt('simulation_input_data_lzk.txt', names=True, dtype=None)

m1 = 1.71079e+06
m2 = 529521
z = 4.6032768
separation = 0.35654176
star_gamma = -0.5160731328052469
vel_disp_1 = 29.998378664129117
vel_disp_2 = 29.998378664129117

#initialize class
evolve_class=EvolveFDFA()

#Calculate timescales
# uses astropy to calculate formation time based on redshift (see below)
# this will output the coalescence_time, not the delay timescale (see paper)
# IMPORTANT NOTE: the redshift at coalescence time needs to be calculated 
# going backwards (like z_at_value in astropy I think)

#age of universe at z=1.75 + delay timescale
coalescence_time = evolve_class(m1, m2, z, separation, star_gamma, vel_disp_1, vel_disp_2)
print("Coalescence Time: ",coalescence_time)

# to get delay timescale, how long it takes from simulation snapshot to coalescence
#we think its sum of below
delay_timescale = evolve_class.t_delay
print("Delay Timescale: ",delay_timescale)

# large scale decay
large_scale_delay=evolve_class.large_scale_decay_time
print("Large Scale Delay: ",large_scale_delay)

# dynamical friction
dynamical_friction=evolve_class.DF_timescale
print("Dynamical Friction: ",dynamical_friction)

# hardening / GW (see the original paper for this construction)
hardening=evolve_class.Hardening_GW_timescale
print("Hardening: ",hardening)

sum_of_LSD_DF_H=large_scale_delay+dynamical_friction+hardening
print("Sum of LSD, DF, & H: ",sum_of_LSD_DF_H)



#plot
#_ = plt.hist2d(np.log10(delay_timescale), np.log10(coalescence_time))
#plt.ylabel(r'log$_{10}(t_{coal}$)', fontsize=16)
#plt.xlabel(r'log$_{10}(t_{delay}$)', fontsize=16)
#plt.show()~
