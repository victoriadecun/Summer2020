import pynbody
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from matplotlib.ticker import NullFormatter
import astropy.units as u

#load simulation
s = pynbody.load('/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096')

h=s.halos()

#convert the units                                                                 
s.physical_units()

def findBH(s):
    BH = s.stars[pynbody.filt.LowPass('tform', 0.0)]
    return BH

#align the galaxy

pynbody.analysis.angmom.faceon(h[5])
BH = findBH(s)
BH_pos = BH['pos']
BHx = BH_pos[:,0]
BHy = BH_pos[:,1]
BHz = BH_pos[:,2]
BH_position = np.array([BHx[0], BHy[0], BHz[0]])

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
x = np.array([vel[0] for vel in velocity])
y = np.array([vel[1] for vel in velocity])
z = np.array([vel[2] for vel in velocity])

#average of these by dividing by total (velocity dispersion)
vel_answer = np.sqrt((x)**2 + (y)**2 + (z)**2)

#divide by total # of stars
velocity = vel_answer.sum() / total_stars

print 'The stars around the black hole in this snapshot are moving at:', velocity
                                                         
