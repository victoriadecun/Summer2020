
import pynbody
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from matplotlib.ticker import NullFormatter

#load snapshot

s=pynbody.load("/media/jillian/h229/h229.cosmo50PLK.3072gst5HbwK1BH.001056/h229.cosmo50PLK.3072gst5HbwK1BH.001056")

#convert the units
s.physical_units()
"""
#  load any available halo
h = s.halos()

# function to find black hole
def findBH(s):
    BHfilter = pynbody.filt.LowPass('tform',0.0)
    BH = s.stars[BHfilter]
    return BH
BH = findBH(s)
print BH

#function to find the halos that the galaxy is in
def findBHhalos(s):
    BH = findBH(s)
    BHhalos = BH['amiga.grp']
    return BHhalos

#using the function the halos
BHhalos = findBHhalos(s)

#printing the halos
print BHhalos

#sorting the halos, indexes/indecis are like an exact address
currenthalo = np.argsort(BHhalos)
print BHhalos[currenthalo]
def getz(s):
    return s.properties['z']

def gettime(s):
    return pynbody.analysis.cosmology.age(s)
for i in currenthalo:

    #which halo are we on?
    currenthalo = BHhalos[i]
    print 'current halo: ', currenthalo

    print i
    
    #put the galaxy you care about in the center of the simulation
    pynbody.analysis.angmom.faceon(h[currenthalo])
    #with pynbody.analysis.halo.center(h[currenthalo], mode='hyb'):

#for halo in BHhalos:
#    halo/galaxy we need to look at is no. 2
#    pynbody.analysis.angmom.faceon(h[halo])
#    prints the black hole
#    print BH

    #this is the position of black hole
    BHposition=BH['pos']

    #printing the position of black hole
    #print BHposition

    #putting the x-values into a column
    BHx= BHposition[[i],0]
    print "x postion", BHx

    #putting the y-values into a column
    BHy= BHposition[[i],1]
    print "y position", BHy

    #putting the z-values into a column
    BHz= BHposition[[i],2]
    print "z position", BHz

    #the .5 is the square root , this is the distance formula
    #distance =((BHx**2)+(BHy**2)+(BHz**2))**(.5)
    #print 'this is the distance :'
    #print "this is the distance :", distance
    data = [currenthalo, BH['iord'][i]] 
    print "This is current halo, ID", data
"""

#  load any available halo                                                                             
h = s.halos()
h4=h[4]

#put galaxy in center of simulation
pynbody.analysis.angmom.faceon(h4)

#create a profile object for stars in 3D                                                               
p = pynbody.analysis.profile.Profile(h[4].s,min=.01,max=2,nbins=50,ndim=3,type='log')


#range out is 0.4, range in 0.04
rout=0.4
rin=0.02
#filter is filt, range of rbins >in and <out
filt=np.where((p['rbins']>rin) & (p['rbins']<rout))
print filt

#show slope of density profile
plt.plot(np.log10(p['rbins']),np.log10(p['density']),'k')
x=np.array(np.log10(p['rbins'][filt]))
y=np.array(np.log10(p['density'][filt]))
m,b=np.polyfit(x,y,1)
plt.plot(x,y,'o')
#plt.xscale('log')
#plt.yscale('log')
#plt.plot.semilogy()
#plt.plot.semilogx()
plt.xlabel('R [kpc]')
plt.ylabel(r'$\rho_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
plt.plot(x,m*x+b)
#plt.twinx()
plt.title('Density Profile of Stars')
plt.show()
print ("This is Slope & Y-Intercept: ", m,b)

#calculate distance for separation

#finds number of BH
def findBH(s):
     BHfilter = pynbody.filt.LowPass('tform',0.0)
     BH = s.stars[BHfilter]
     return BH
BH = findBH(s)
print('number of black holes: ',len(BH))

#how far away BH from galaxy center
#it's x,y,z positions
with pynbody.analysis.halo.center(h[4], mode='hyb'):
    print (h[4]['pos'][0])
    print (h[4]['pos'][1])
    print (h[4]['pos'][2])

print('pos')

#position of BH
BHposition=BH['pos']
print('this is BH position: ', BHposition)

#put x-values in column
BHx=BHposition[:,0]
print('this is the x-position: ', BHx)

#put y-values in column
BHy=BHposition[:,1]
print('this is the y-position: ', BHy)

#put z-values in column
BHz=BHposition[:,2]
print('this is the z-position: ', BHz)

#.5 is the square root; distance
distance=((BHx**2)+(BHy**2)+(BHz**2))**(.5)
print('this is the distance: ', distance)

#Calculate velocity dispersions

BH_position=np.array([BHx[0],BHy[0],BHz[0]])
#this puts sphere around BH - 0.116 is how i got 99 stars 
radius= 0.116  #kpc
sphere = pynbody.filt.Sphere(radius, cen =(BH_position))
stars = s.stars[0:]
in_sphere = stars[sphere]
total_stars = len(in_sphere)
print("Total stars: ",total_stars)

#find velocity
velocity=in_sphere['vel']

#find velocity of stars in x,y,z
x = np.std(velocity[0])
y = np.std(velocity[1])
z = np.std(velocity[2])

#average of these by dividing by total (velocity dispersion)
vel_answer = np.sqrt((x)**2 + (y)**2 + (z)**2)

print("Velocity Dispersion: ",vel_answer)

#Calculating delay time

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)

from scipy.interpolate import interp1d, interp2d
from scipy.special import gamma as Gamma_Function
from scipy.integrate import quad
from scipy.special import hyp2f1
from scipy import constants as ct
import numpy as np
import pdb

from mbhbinaries import (
    mass_ratio_func,
    MassiveBlackHoleBinaries,
    AnalyticApproximations,
)

Msun = 1.989e30


class EvolveFDFA(MassiveBlackHoleBinaries):
    # put all quantities into arrays to determine which is major galaxy and minor galaxy'

    def __init__(self, e_0=0.0):
        """
        data = np.genfromtxt(fname, names=True, dtype=None)


        evolve_key_dict = {
            "m1": "mass_new_prev_in",
            "m2": "mass_new_prev_out",
            "z": "redshift",
            "separation": "separation",
            "star_gamma": "star_gamma",
            "vel_disp_1": "vel_disp_prev_in",
            "vel_disp_2": "vel_disp_prev_out",
        }

        for key, col_name in evolve_key_dict.items():
            setattr(self, key, data[col_name])
        """

    def __call__(self, m1, m2, z, separation, star_gamma, vel_disp_1, vel_disp_2):

        self.m1, self.m2 = 4.93481e+06, 3.36428e+06
        self.z = 1.7536451
        self.separation = 3.38493483e+02
        self.star_gamma = -0.5957736982667012
        self.vel_disp_1, self.vel_disp_2 = 189.03323719, vel_disp_2

        # find index of major and minor
        major_1 = self.m1 >= self.m2
        major_2 = self.m1 < self.m2

        # major black hole mass
        self.M = self.m1 * major_1 + self.m2 * major_2
        # minor black hole mass
        self.m = self.m1 * major_2 + self.m2 * major_1

        self.m1 = self.M
        self.m2 = self.m

        # small s denotes secondary,small m is primary (same as paper)
        self.vel_disp_m = self.vel_disp_1 * major_1 + self.vel_disp_2 * major_2
        self.vel_disp_s = self.vel_disp_1 * major_2 + self.vel_disp_2 * major_1

        # find large scale orbital decay time
        self.gamma = np.clip(self.star_gamma, 0.55, 2.49)

        # r_infl determined analytically
        self.r_infl = AnalyticApproximations.influence_radius(self.M, self.vel_disp_m)

        self.R_e_m = self.separation

        self.q = mass_ratio_func(self.M, self.m)

        self.e_0 = 0.0  # e_0

        self.Lambda = np.clip(
            2 ** (3.0 / 2.0) * (self.vel_disp_m / self.vel_disp_s),
            np.exp(2.0),
            np.exp(6.0),
        )

        self.b = self.gamma - 3.0 / 2.0

        self.formation_time = cosmo.age(self.z).value * 1e9
        self.t_delay, self.e_f = self.calculate_timescale()
        self.coalescence_time = self.formation_time + self.t_delay

        return self.coalescence_time
        print ("Coalescence Time: ", self.coalescence_time)

    def calculate_timescale(self, return_arr=False):
        self.FD_FA_large_scale_orbital_decay_timescale()
        self.FD_FA_Dynamical_Friction_timescale()
        self.FD_FA_hardening_timescale()

        # TODO: add eccentricity capabilities
        self.e_f = self.e_0

        if return_arr:
            return np.asarray(
                [
                    self.large_scale_decay_time,
                    self.DF_timescale,
                    self.Hardening_GW_timescale,
                ]
            ).T

        return (
            self.large_scale_decay_time
            + self.DF_timescale
            + self.Hardening_GW_timescale,
            self.e_f,
        )
        print ("Timescale: ", self.large_scale_decay_time + self.DF_timescale + self.Hardening_GW_timescale, self.e_f)

    def T_gx_star_1(self):
        return 1e9 * (
            0.06
            * (2.0 / np.log(self.Lambda))
            * (self.R_e_m / 10.0) ** 2
            * (self.vel_disp_m / 300.0)
            * (1e8 / self.m)
        )  # yr

    def T_gx_star_2(self):
        return 1e9 * (
            0.15
            * (2.0 / np.log(self.Lambda))
            * (self.R_e_m / 10.0)
            * (self.vel_disp_m / 300.0) ** 2
            * (100.0 / self.vel_disp_s) ** 3
        )  # yr

    def FD_FA_large_scale_orbital_decay_timescale(self):
        # Equations 54, 56, 57
        out = np.asarray([self.T_gx_star_1(), self.T_gx_star_2()]).T
        self.large_scale_decay_time = np.max(out, axis=1)
        return

    def alpha_func(self):
        self.alpha = (
            (Gamma_Function(self.gamma + 1) / Gamma_Function(self.gamma - 1 / 2))
            * (4.0 / 3.0)
            * np.pi ** (-1 / 2)
            * 2.0 ** (self.b - self.gamma)
            * self.ksi ** 3
            * hyp2f1(3.0 / 2.0, -self.b, 5.0 / 2.0, self.ksi ** 2 / 2.0)
        )
        return

    def beta_func(self):
        beta_integral_vectorized = np.frompyfunc(beta_integral_func, 2, 1)
        integral = beta_integral_vectorized(self.ksi, self.b).astype(np.float64)

        self.beta = (
            (Gamma_Function(self.gamma + 1) / Gamma_Function(self.gamma - 1 / 2))
            * 4
            * np.pi ** (-1 / 2)
            * 2 ** -self.gamma
            * integral
        )
        return

    def delta_func(self):
        self.delta = (
            (Gamma_Function(self.gamma + 1) / Gamma_Function(self.gamma - 1 / 2))
            * 8
            * np.pi ** (-1 / 2)
            * (2 ** (-self.gamma - 1) / (self.b + 1))
            * self.ksi
            * (0.04 ** (self.b + 1) - (2 - self.ksi ** 2) ** (self.b + 1))
        )
        return

    def T_dot_bare(self):
        # in T_dot_bare the coulomb logarith is set to 6.0
        return (
            1.5e7
            * (6.0 * self.alpha + self.beta + self.delta) ** -1
            / ((3 / 2.0 - self.gamma) * (3.0 - self.gamma))
            * (self.chi ** (self.gamma - 3.0 / 2.0) - 1)
            * (self.M / 3e9) ** (1 / 2)
            * (self.m / 1e8) ** -1
            * (self.r_infl / 300) ** (3 / 2)
        )

    def T_dot_gx(self):
        return (
            1.2e7
            * (np.log(self.Lambda) * self.alpha + self.beta + self.delta) ** -1
            / ((3.0 - self.gamma) ** 2)
            * (self.chi ** (self.gamma - 3.0) - 1)
            * (self.M / 3e9)
            * (100 / self.vel_disp_s) ** 3
        )  # years

    def find_a_crit(self):
        return self.r_infl * (self.m / (2 * self.M)) ** (1 / (3 - self.gamma))  # pc

    def FD_FA_Dynamical_Friction_timescale(self, use_interp=False):
        # a_crit = find_a_crit(r_infl, m, M, gamma)
        # chi = a_crit/r_infl
        # self.find_a_h()
        # self.chi = self.a_h/self.r_infl

        self.find_chi()

        self.ksi = 1.0
        self.alpha_func()
        self.beta_func()
        self.delta_func()

        out = np.asarray([self.T_dot_bare(), self.T_dot_gx()]).T

        self.DF_timescale = np.min(out, axis=1)
        return

    """
    # Eccentricity related functions
    def k_func(M, m):
        return 0.6 + 0.1*np.log10((M + m) / 3e9)

    def p_e_func(e, M, m):
        k = k_func(M,m)
        return (1-e**2)*(k + (1-k) * (1-e**2)**4)
    """

    def FD_FA_hardening_timescale(self, psi=0.3, phi=0.4):
        #  Equations 61-63
        # includes gravitational regime

        # Below timescale is from FD's paper
        # if e ==0:
        #    p_e = 1.0
        # else:
        # p_e = p_e_func(e, M, m)
        # T_h_GW = (1.2e9 * (r_infl/300.)**((10 + 4*psi)/(5 + psi))
        #   * ((M+m)/3e9)**((-5-3*psi)/(5+psi)) * phi**(-4/(5+psi))
        #   * (4*q/(1+q)**2)**((3*psi - 1)/(5 + psi))* p_e) #years

        # We decided to use Vasiliev's eq. 25

        # f_e = self.f_e_func() * (self.e_f !=0.0) + 1.0 * (self.e_f==0.0)
        f_e = 1.0

        T_h_GW = (
            1.7e8
            * (self.r_infl / 30.0) ** ((10 + 4 * psi) / (5 + psi))
            * ((self.M + self.m) / 1e8) ** ((-5 - 3 * psi) / (5 + psi))
            * phi ** (-4 / (5 + psi))
            * (4 * self.q / (1 + self.q) ** 2) ** ((3 * psi - 1) / (5 + psi))
            * f_e ** ((1 + psi) / (5 + psi))
            * 20 ** psi
        )  # years

        self.T_GW = self.find_T_GW()
        self.Hardening_GW_timescale = T_h_GW * (self.q >= 1e-3) + self.T_GW * (
            self.q < 1e-3
        )
        return

    def find_chi(self):
        # combines equation 3b from Merritt et al 2009 and hardening timescale from merritt 2013
        self.chi = 0.25 * (self.q / (1.0 + self.q))
        return

    def find_a_h(self):
        self.a_h = (
            36.0
            * (self.q / (1 + self.q) ** 2)
            * ((self.M + self.m) / 3e9)
            * (self.vel_disp_m / 300) ** -2
        )  # pc
        return

    """
    # Eccentricity related
    def f_e_func(self):
        return (1-self.e_f**2)**(7/2)/(1 + (73./24.)*self.e_f**2 + (37.0/96.)*self.e_f**4)
    """

    def find_T_GW(self):
        self.find_a_h()
        m1 = self.M * Msun * ct.G / ct.c ** 2
        m2 = self.m * Msun * ct.G / ct.c ** 2
        beta = 64.0 / 5.0 * m1 * m2 * (m1 + m2)
        T_GW_meters = self.a_h ** 4 / beta
        T_GW_seconds = T_GW_meters / ct.c
        T_GW_yr = T_GW_seconds / (ct.Julian_year)
        return T_GW_yr


def beta_integrand_func(x, ksi, b):
    return x ** 2 * (2 - x ** 2) ** b * np.log((x + ksi) / (x - ksi))


def beta_integral_func(ksi, b):
    integral, err = quad(beta_integrand_func, ksi, 1.4, args=(ksi, b))
    return integral
