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
data = [5, BH['iord'], gettime(s),getz(s), BH['mass'], BH['r'], starmass, gasmass, virialmass]

print('this is the data: ', data)

