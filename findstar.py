import pynbody
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from matplotlib.ticker import NullFormatter

#find how many stars & load 

s=pynbody.load("/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096")

#convert the units
s.physical_units()


#this finds stars
#stars=s.stars[0:]
#print(stars)

#how many stars in galaxy 5 (halo 5)
h=s.halos()
h5=h[5]

#put galaxy in center of simulation
pynbody.analysis.angmom.faceon(h5)

#finds gas, dark matter and stars
#stars=s.stars[0:]
#print('ngas = %e, ndark = %e, nstar = %e\n'%(len(h5.gas),len(h5.dark),len(h5.star)))

#finds number of BH
def findBH(s):
     BHfilter = pynbody.filt.LowPass('tform',0.0)
     BH = s.stars[BHfilter]
     return BH
BH = findBH(s)
print('number of black holes: ',len(BH))

#how far away BH from galaxy center
#it's x,y,z positions
with pynbody.analysis.halo.center(h[5], mode='hyb'):
    print (h[5]['pos'][0])
    print (h[5]['pos'][1])
    print (h[5]['pos'][2])

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

#find mass and ID number
def getz(s):
    return s.properties['z']

def gettime(s):
    return pynbody.analysis.cosmology.age(s)

starmass = h[5].s['mass'].sum()
gasmass = h[5].g['mass'].sum()
virialmass = starmass+gasmass+h[5].d['mass'].sum()
data = ['halo: ',5,'BH_ID: ', BH['iord'],'time: ', gettime(s),'z: ',getz(s),'BH_mass: ', BH['mass'],'Radius: ', BH['r'],'star_mass: ', starmass,'gas_mass: ', gasmass,'virial_mass: ', virialmass]

print('this is the data: ', data)

"""
#make an image of galaxy
pynbody.plot.stars.render(s,width='10 kpc')
plt.show()
"""

#density profiles
#create a profile object for stars in 3D
#p = pynbody.analysis.profile.Profile(h[5].s,min=.01,max=2,nbins=50,ndim=3,type='log')

# make a 3D density plot of the dark matter (note ndim=3 in the constructor below)                                                             
#p1 = pynbody.analysis.profile.Profile(h[5].d,min=.01,max=2,nbins=50,ndim=3,type='log')

#make a 3D density plot of gas                                                                                                                 
#p2 = pynbody.analysis.profile.Profile(h[5].g,min=.01,max=2,nbins=50,ndim=3,type='log')

#make a 3D density plot of all combined                                                                                                        
#p3 = pynbody.analysis.profile.Profile(h[5],min=.01,max=2,nbins=50,ndim=3,type='log')

# make the figure and sub plots                                                                                                                
#f, axs = plt.subplots(2,2,figsize=(20,6))

# make the plot

"""
axs[0,0].plot(p['rbins'],p['density'], 'k')
axs[0,0].semilogy()
axs[0,0].semilogx()
axs[0,0].set_xlabel('R [kpc]')
axs[0,0].set_ylabel(r'$\rho_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
axs[0,0].set_title('Density Profile of Stars')

axs[1,0].plot(p1['rbins'],p1['density'], 'k')
axs[1,0].semilogy()
axs[1,0].semilogx()
axs[1,0].set_xlabel('R [kpc]')
axs[1,0].set_ylabel(r'$\rho_{DM}$ [M$_{\odot}$ kpc$^{-3}$]')
axs[1,0].set_title('Density Profile of Dark Matter')

axs[0,1].plot(p2['rbins'],p2['density'], 'k')
axs[0,1].semilogy()
axs[0,1].semilogx()
axs[0,1].set_xlabel('R [kpc]')
axs[0,1].set_ylabel(r'$\rho_{G}$ [M$_{\odot}$ kpc$^{-3}$]')
axs[0,1].set_title('Density Profile of Gas')

axs[1,1].plot(p3['rbins'],p3['density'], 'k')
axs[1,1].semilogy()
axs[1,1].semilogx()
axs[1,1].set_xlabel('R [kpc]')
axs[1,1].set_ylabel(r'$\rho_{ALL}$ [M$_{\odot}$ kpc$^{-3}$]')
axs[1,1].set_title('Density Profile of Stars, Dark Matter, & Gas')

plt.subplots_adjust(left=None,bottom=None,right=None,top=0.9,wspace=0.5,hspace=0.5)
plt.show()
"""

#make slope of stellar density profile only

#create a profile object for stars in 3D                                                                
p = pynbody.analysis.profile.Profile(h[5].s,min=.1,max=1,nbins=50,ndim=3,type='log')

plt.plot(p['rbins'],p['density'],'k')
x=np.array(p['rbins'])
y=np.array(p['density'])
m,b=np.polyfit(x,y,1)
plt.plot(x,y,'o')
plt.xscale('log')
plt.yscale('log')
#plt.plot.semilogy()
#plt.plot.semilogx()
plt.xlabel('R [kpc]')
plt.ylabel(r'$\rho_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
plt.plot(x,m*x+b)
plt.twinx()
plt.title('Density Profile of Stars')
plt.show()
